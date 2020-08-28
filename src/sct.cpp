#include <vector>
#include <math.h>
#include "titanlib.h"
#include <assert.h>
#include <iostream>
#include <numeric>
#include <exception>
#include <gsl/gsl_multimin.h>
#include <gsl/gsl_linalg.h>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/vector_proxy.hpp>
#include <boost/numeric/ublas/triangular.hpp>
#include <boost/numeric/ublas/lu.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/ublas/matrix.hpp>


// helpers
bool invert_matrix (const boost::numeric::ublas::matrix<float>& input, boost::numeric::ublas::matrix<float>& inverse);
ivec remove_flagged(ivec indices, ivec flags);
vec compute_vertical_profile(const vec& elevs, const vec& oelevs, const vec& values, int num_min_prof, double min_elev_diff);
double basic_vertical_profile_optimizer_function(const gsl_vector *v, void *data);
vec basic_vertical_profile(const int n, const double *elevs, const double t0, const double gamma);
double vertical_profile_optimizer_function(const gsl_vector *v, void *data);
vec vertical_profile(const int n, const double *elevs, const double t0, const double gamma, const double a, const double h0, const double h1i);

// start SCT //
ivec titanlib::sct(const vec& lats,
        const vec& lons,
        const vec& elevs,
        const vec& values,
        const ivec& obs_to_check,
        const vec& background_values,
        std::string background_elab_type,
        int num_min,
        int num_max,
        float inner_radius,
        float outer_radius,
        int num_iterations,
        int num_min_prof,
        float min_elev_diff,
        float min_horizontal_scale,
        float max_horizontal_scale,
        int kth_closest_obs_horizontal_scale,
        float vertical_scale,
        float value_min,
        float value_max,
        float sig2o_min,
        float sig2o_max,
        const vec& eps2,
        const vec& tpos_score,
        const vec& tneg_score,
        const vec& t_sod,
        bool debug,
        vec& score,
        vec& rep,
        vec& sod,
        vec& num_inner,
        vec& horizontal_scale,
        vec& an_inc,
        vec& an_res,
        vec& cv_res,
        vec& innov,
        vec& idi,
        vec& idiv,
        vec& sig2o) {
/*

 -- Spatial Consistency Test --

 Description: 
  Flag observations that are (likely) affected by gross measurement errors (GE) based on neighbouring observations. Observations are affected by GEs if their values are (a) [not related to the actual atmospheric state] OR (b) [affected by such large representativeness errors (REs) that they are difficult to reconstruct using the neighbouring observations].

 Reference: 
  Lussana, C., Uboldi, F. and Salvati, M.R. (2010), A spatial consistency test for surface observations from mesoscale meteorological networks. Q.J.R. Meteorol. Soc., 136: 1075-1088. doi:10.1002/qj.622 (Reference Abbreviation: LUS10)
 See also the wiki-pages https://github.com/metno/titanlib/wiki/Spatial-consistency-test

 Definitions:
 + obs_to_check(i) = 1 check ith observation, otherwise use it but do not check it
 + Observations (p_outer-vector): 
    (observed_)values = true_values + observation_errors
    NOTE: true_values are unknown, useful for theory
 + Background (p_outer-vector):
    background_value = true_value + background_errors
 + Observation_errors (p_outer-vector): 
   random variable that follows a multivariate normal (MVN) distribution with:
     R = p_outer x p_outer covariance matrix, R = sig2o * identity_matrix 
     mean (p_outer-vector) = 0 (for all s elements)
 + Background_errors (p_outer-vector):
   random variable that follows a MVN distribution with:
     S= p_outer x p_outer covariance matrix, S = exponential 
     mean (p_outer-vector) = 0
 + Innovation (p_outer-vector)= observations - background 
 + Analysis residuals (p_outer-vector)= observations - analysis 
 + Analysis increments (p_outer-vector)= analysis - background
 + Cross-validation indicates leave-one-out cross-validation
 + (SCT) Score (score, p_outer-vector)
   score = cv-analysis residuals * analysis residuals / sig2o, LUS10 Eq.(22)
 + Coefficient of representativeness (corep, p_outer-vector)
   rep = innovation * analysis residuals / sig2o, LUS10 Eq.(32) numerator
 + Spatial Outlier Detection (SOD) score (sod, p_outer-vector)
    used to distinguish between acceptable REs and GEs
     sod = (deviation from SCT-score areal average)^2 / spread 

 Algorithm: 
  Loop over the s observations (observation index is "curr"). Select the subset of p_outer closest observations to the curr-th observation, then perform SCT considering only this subset.
  curr-th observation is flagged if curr-th score is a) larger than a pre-set threshold AND b) an outlier compared to the statistics of score inside the inner circle
  condition b) outlier if sod>threshold. Helps to avoid flagging REs as GEs, at least in those regions where enough observations are avaialable

 Returned values:
  flags. -9999 = not checked; 0 = passed (good); 1 = failed (bad); 11 = isolated (<2 inside inner); 12 = isolated (<num_min inside outer)

*/
    const int p = values.size();
    if( lats.size() != p || lons.size() != p || elevs.size() != p || values.size() != p || tpos_score.size() != p || tneg_score.size() != p || eps2.size() != p)
        throw std::runtime_error("Dimension mismatch");
    if(num_min < 2)
        throw std::invalid_argument("num_min must be > 1");
    if(num_max < num_min)
        throw std::invalid_argument("num_max must be > num_min");
    if(num_iterations < 1)
        throw std::invalid_argument("num_iterations must be >= 1");
    if(min_elev_diff <= 0)
        throw std::invalid_argument("min_elev_diff must be > 0");
    if(min_horizontal_scale <= 0)
        throw std::invalid_argument("min_horizontal_scale must be > 0");
    if(kth_closest_obs_horizontal_scale <= 0)
        throw std::invalid_argument("kth_closest_obs_horizontal_scale must be > 0");
    if(vertical_scale <= 0)
        throw std::invalid_argument("vertical_scale must be > 0");
    if(inner_radius < 0)
        throw std::invalid_argument("inner_radius must be >= 0");
    if(outer_radius < inner_radius)
        throw std::invalid_argument("outer_radius must be >= inner_radius");
    if(value_max <= value_min)
        throw std::invalid_argument("value_max must be > value_min");
    if(sig2o_max <= sig2o_min)
        throw std::invalid_argument("sig2o_max must be > sig2o_min");
    if(sig2o_min < 0)
        throw std::invalid_argument("sig2o_min must be >= 0");
    if(background_elab_type != "vertical_profile" && background_elab_type != "mean_outer_circle" && background_elab_type != "external")
        throw std::invalid_argument("background_elab_type must be one of vertical_profile, mean_outer_circle or external");
    if(background_elab_type == "vertical_profile" && num_min_prof<0)
        throw std::invalid_argument("num_min_prof must be >=0");
    if( (obs_to_check.size() != p && obs_to_check.size() != 1 && obs_to_check.size() !=0) ) 
        throw std::invalid_argument("'obs_to_check' has an invalid length");
    if(background_elab_type == "external" &&  background_values.size() != p)
        throw std::runtime_error("Background vector dimension mismatch");

    // initializations
    float na = -9999.; // code for "Not Available". Any type of missing data
    ivec flags(p, na);
    score.clear();
    score.resize(p, na);
    sod.clear();
    sod.resize(p, na);
    num_inner.clear();
    num_inner.resize(p, na);
    horizontal_scale.clear();
    horizontal_scale.resize(p, na);
    an_inc.clear();
    an_inc.resize(p, na);
    an_res.clear();
    an_res.resize(p, na);
    cv_res.clear();
    cv_res.resize(p, na);
    innov.clear();
    innov.resize(p, na);
    idi.clear();
    idi.resize(p, na);
    idiv.clear();
    idiv.resize(p, na);
    sig2o.clear();
    sig2o.resize(p, na);
    rep.clear();
    rep.resize(p, na);


    // preliminary range check
    int thrown_out = 0; 
    for(int j=0; j < p; j++) {
        if((values[j] < value_min) || (values[j] > value_max)) { 
          flags[j] = 1;
          thrown_out++;
        }
    }
    std::cout << "Range check is removing " << thrown_out << " observations " << std::endl;

    // if obs_to_check is empty then check all
    bool check_all = obs_to_check.size() != p;

    // KDtree has to do with fast computation of distances
    titanlib::KDTree tree(lats, lons);

    // SCT iterations
    for(int iteration = 0; iteration < num_iterations; iteration++) {
        double s_time0 = titanlib::util::clock();
        
        // reset this number each loop (this is for breaking if we don't throw anything new out)
        int thrown_out = 0; 
        
        // diagnostic. count the number of times OI is performed
        int count_oi = 0;

        // loop over observations
        for(int curr=0; curr < p; curr++) {

            if(debug) std::cout << "===> curr " << curr << "===============" << std::endl;

            // observation not to be checked
            if(!check_all && obs_to_check[curr] != 1) {
                if(debug) std::cout << "..not to check " << curr << std::endl;
                continue;
            }

            // break out if station already flagged
            if( flags[curr] == 1 || flags[curr] == 0) {
                if(debug) std::cout << "..checked " << curr << std::endl;
                continue;
            }

            // if no neighbours inside the inner circle, then distinction between RE and GE not reliable
            vec distances;
            ivec neighbour_indices = tree.get_neighbours_with_distance(lats[curr], lons[curr], inner_radius, 2, true, distances);
            neighbour_indices = remove_flagged(neighbour_indices, flags);
            if(neighbour_indices.size() < 2) {
                flags[curr] = 11;
                if(debug) std::cout << "@@isolated (inner) " << curr << std::endl;
                continue;
            }
 
            // get all neighbours that are close enough (inside outer circle)
            neighbour_indices = tree.get_neighbours_with_distance(lats[curr], lons[curr], outer_radius, num_max, true, distances);
            neighbour_indices = remove_flagged(neighbour_indices, flags);
            if(neighbour_indices.size() < num_min) {
                flags[curr] = 12;
                if(debug) std::cout << "@@isolated (outer) " << curr << std::endl;
                continue; 
            }
            int p_inner = neighbour_indices.size();

            // call SCT with this box(=outer circle) 
            vec lons_box = titanlib::util::subset(lons, neighbour_indices);
            vec elevs_box = titanlib::util::subset(elevs, neighbour_indices);
            vec lats_box = titanlib::util::subset(lats, neighbour_indices);
            vec values_box = titanlib::util::subset(values, neighbour_indices);
            vec eps2_box = titanlib::util::subset(eps2, neighbour_indices);
            int p_outer = neighbour_indices.size();
            if(debug) std::cout << "p_outer p_inner " << p_outer << " " << p_inner << " " << std::endl;

            // Compute the background, if needed
            vec bvalues_box;
            if( background_elab_type == "vertical_profile") {
                bvalues_box = compute_vertical_profile(elevs_box, elevs_box, values_box, num_min_prof, min_elev_diff);
            } else if( background_elab_type == "mean_outer_circle"){
                double mean_val = std::accumulate(values_box.begin(), values_box.end(), 0.0) / values_box.size();
                bvalues_box.resize(p_outer, na);
                for(int i = 0; i < p_outer; i++) {
                    bvalues_box[i] = mean_val;
                }
            } else if( background_elab_type == "external"){
                bvalues_box = titanlib::util::subset(background_values, neighbour_indices);
            } 
            if(debug) std::cout << "background ok" << std::endl;
            
            // Compute innovations (= observations - background)
            boost::numeric::ublas::vector<float> d(p_outer);
            bool small_innov = true;
            for(int i=0; i<p_outer; i++) {
                // check the background is in a range of acceptable values
                if(bvalues_box[i] < value_min) bvalues_box[i] = value_min;
                if(bvalues_box[i] > value_max) bvalues_box[i] = value_max;
                // innovations
                d(i) = values_box[i] - bvalues_box[i];
                if ( d(i) > sig2o_min) small_innov = false;
            }
            
            /* if observations and background are almost identical, then take a shortcut
               flag = 0 for all observations in the inner circle
               NOTE: when innovations are all exectly =0, then the routine crashes
             */
            if (small_innov) {
                for(int i = 0; i < p_outer; i++) {
                    int index = neighbour_indices[i];
                    if ( flags[index] == na) {
                        float dist = distances[i];
                        if(dist <= inner_radius) {
                            score[index] = 0;
                            rep[index] = 0;
                            sod[index] = 0;
                            num_inner[index] = p_inner;
                            horizontal_scale[index] = na;
                            an_inc[index] = na;
                            an_res[index] = na;
                            cv_res[index] = na;
                            innov[index] = d(i);
                            idi[index] = na;
                            idiv[index] = na;
                            sig2o[index] = na;
                            flags[index] = 0;
                            if(debug) {
                                int j = index;
                                std::cout << std::setprecision(3) << " small_innov - flag=0 - index innov ares ainc cvres idi idiv: " << j << " " << innov[j] << " " << an_res[j] << " " << an_inc[j] << " " << cv_res[j] << " " << idi[j] << " " << idiv[j] << " "  << std::endl;
                            }
                        }
                    }
                }
                continue; // jump to the next observation
            }
 
            /* Compute Dh. 
             The location-dependent horizontal de-correlation lenght scale 
             used for the background error correlation matrix */
            boost::numeric::ublas::matrix<float> disth(p_outer, p_outer);
            boost::numeric::ublas::matrix<float> distz(p_outer, p_outer);
            boost::numeric::ublas::vector<float> Dhbox(p_outer);

            for(int i=0; i < p_outer; i++) {
                vec Dhbox_vector(p_outer);
                for(int j=0; j < p_outer; j++) {
                    disth(i, j) = titanlib::util::calc_distance(lats_box[i], lons_box[i], lats_box[j], lons_box[j]);
                    distz(i, j) = fabs(elevs_box[i] - elevs_box[j]);
                    if(i != j) {
                        if(i < j)
                            Dhbox_vector[j - 1] = disth(i, j);
                        else if(i > j)
                            Dhbox_vector[j] = disth(i, j);
                    }
                }
                // Dhbox(i) = titanlib::util::compute_quantile(0.10, Dhbox_vector);  // quantile
                Dhbox(i) = titanlib::util::findKclosest( kth_closest_obs_horizontal_scale, Dhbox_vector); // k-th closest observations
            }

            double Dhbox_mean = std::accumulate(std::begin(Dhbox), std::end(Dhbox), 0.0) / Dhbox.size();
            if(Dhbox_mean < min_horizontal_scale) {
                if(debug) std::cout << "Dhbox_mean (<min_horizontal_scale) " << Dhbox_mean << std::endl;
                Dhbox_mean = min_horizontal_scale;
            }
            if(Dhbox_mean > max_horizontal_scale) {
                if(debug) std::cout << "Dhbox_mean (>max_horizontal_scale) " << Dhbox_mean << std::endl;
                Dhbox_mean = max_horizontal_scale;
            }
            if(debug) std::cout << "Dhbox_mean " << Dhbox_mean << std::endl;
            
            // Compute S + eps2*I and store it in S 
            boost::numeric::ublas::matrix<float> S(p_outer,p_outer);
            boost::numeric::ublas::matrix<float> Sinv(p_outer,p_outer);
            for(int i=0; i < p_outer; i++) {
                for(int j=0; j < p_outer; j++) {
                    double value = std::exp(-.5 * std::pow((disth(i, j) / Dhbox_mean), 2) - .5 * std::pow((distz(i, j) / vertical_scale), 2));
                    if(i==j) { // weight the diagonal, this also ensure an invertible matrix
                        value = value + eps2_box[i];
                    }
                    S(i,j) = value;
                }
            }
            
            /* ---------------------------------------------------
               Beginning of real SCT
               ------------------------------------------------------*/
            // Compute ( S + eps2 * I )^(-1) ans store it in Sinv
            bool b = invert_matrix(S, Sinv);
            if(b != true) {
                // TODO: flag differently or give an error???
                // NOTE: should never happen, since S + eps2*I is made of columns that are all linearly independent
                //       (this is also why is convenient to use this formulation)
                if(debug) std::cout << "ooo Problem in matrix inversion ooo " << curr << std::endl;
                continue;
            }

            // Definitions
            boost::numeric::ublas::vector<float> Zinv(p_outer), Sinv_d(p_outer), Sinv_1(p_outer), ainc(p_outer), ares(p_outer), cvres(p_outer), rep_temp(p_outer), idibox(p_outer), idivbox(p_outer);

            // Matrix multiplications
            for(int i=0; i<p_outer; i++) {
                // Recover the actual S from S + eps2*I (unweight the diagonal)
                S(i,i) -= eps2_box[i];
                // Sinv_d = ( S + eps2 * I )^(-1) * innovation
                // Sinv_1 = row sums of ( S + eps2 * I )^(-1) (used for IDI)
                double acc = 0;
                double acc1 = 0;
                for(int j=0; j<p_outer; j++) {
                    acc += Sinv(i,j)*d(j);
                    acc1 += Sinv(i,j);
                }
                Sinv_d(i) = acc;
                Sinv_1(i) = acc1;
            }
           
            //    analysis increment = analysis - background = S * ( S + eps2 * I )^(-1) * innovation
            //                  Zinv = inverse of the elements on the diagonal of ( S + eps2 * I )^(-1)
            //    analysis residuals = observations - analysis (= innovation - analysis increment)
            // cv-analysis residuals = observations - cv-analysis, LUS10 Eq.(A13)
            // observ error variance = mean( analysis residuals * innovation), LUS10 Eq.(32)
            //                   sod = spatial outlier detection score
            boost::numeric::ublas::vector<float> scorebox_numerator(p_outer);
            double M2a = 0;
            double M2b = 0;
            double scorebox_mean = 0;      // scorebox mean
            double sig2obox = 0;           // sig2obox mean
            for(int i=0; i<p_outer; i++) {
                double acc = 0;
                double acc1 = 0;
                for(int j=0; j<p_outer; j++) {
                    acc += S(i,j)*Sinv_d(j);
                    acc1 += S(i,j)*Sinv_1(j);
                }
                idibox(i) = acc1;
                ainc(i) = acc;
                Zinv(i) = 1 / Sinv(i,i); 
                ares(i) = d(i) - ainc(i);
                cvres(i) = Zinv(i) * Sinv_d(i); 
                idivbox(i) = 1 - Zinv(i) * Sinv_1(i); 
                // (Welford's online algorithm for mean and variance)
                // estimation of average score inside outer circle
                scorebox_numerator(i) = cvres(i) * ares(i);
                float delta = scorebox_numerator(i) - scorebox_mean;
                scorebox_mean += delta / (i+1);
                float delta2 = scorebox_numerator(i) - scorebox_mean;
                M2a += delta * delta2;
                // estimation of sig2o (observation error variance) 
                rep_temp(i) = fabs( d(i) * ares(i)); // I've added fabs(...) 
                delta = rep_temp(i) - sig2obox;
                sig2obox += delta / (i+1);
                delta2 = rep_temp(i) - sig2obox;
                M2b += delta * delta2;
            }

            // finalize estimation of sig2o sample variance
            double sig2obox_var  = M2b / (p_outer-1); 
            // sig2o variance of mean (see e.g. Taylor "an Intro to error analysis..." 1982, p. 102)
            double sig2obox_mean_var = sig2obox_var / p_outer;
            // one should have an idea of the range of allowed values for sig2o
            if(sig2obox < sig2o_min) sig2obox = sig2o_min;
            if(sig2obox > sig2o_max) {
              sig2obox = sig2o_max;
              sig2obox_mean_var = 0; // high enough uncertainty is included in sig2o_max
            }
            if(debug) std::cout << std::setprecision(3) << "sig2obox (var) " << sig2obox << " " << sig2obox_mean_var << " " << std::endl;
 
            // finalize scorebox_mean
            scorebox_mean = scorebox_mean / sig2obox;
            // scorebox sample variance
            double scorebox_var  = M2a / (p_outer-1) * 1/(sig2obox*sig2obox);
            // scorebox variance of mean (see e.g. Taylor "an Intro to error analysis..." 1982, p. 102)
            double scorebox_mean_var = scorebox_var / p_outer;
            if(debug) std::cout << std::setprecision(3) << "scorebox_mean (var) " << scorebox_mean << " " << scorebox_mean_var << " " << std::endl;

            /* score = cv-analysis residuals * analysis residuals / sig2o, LUS10 Eq.(22)
               rep = innovation * analysis residuals / sig2o, LUS10 Eq.(32) numerator
               i-th observation is flagged if i-th score is a) larger than a pre-set threshold AND b) an outlier compared to the statistics of score inside the inner circle
               condition b) helps to avoid flagging REs as GEs, at least in those regions where enough observations are avaialable
               sod (spatial outlier detection score) = (deviation from the mean)^2 / (spread + uncertainty in the mean)
               outlier if sod>threshold 
               */

            // update flags inside the inner circle

            // step 1. select observations possibly to flag and determine max score in the box
            int index_scorebox_max = -1;
            int i_scorebox_max = -1;
            float scorebox_max = -1;
            int count_updates = 0;
            ivec update(p_outer, 0);  // Keep track of which observations need an update on the flags
            boost::numeric::ublas::vector<float> scorebox(p_outer), sodbox(p_outer), corep(p_outer);
            for(int i = 0; i < p_outer; i++) {
                int index = neighbour_indices[i];
                if ( flags[index] == na) {
                    float dist = distances[i];
                    if(dist <= inner_radius) {
                        scorebox(i) = scorebox_numerator(i) / ( sig2obox + sig2obox_mean_var);
                        sodbox(i) = std::pow( ( scorebox_numerator(i) - scorebox_mean), 2) / (scorebox_var + scorebox_mean_var);
                        corep(i) = rep_temp(i) / sig2obox;
                        // condition defining bad observations 
                        if( ((cvres(i) > 0 && scorebox(i) > tpos_score[index]) || (cvres(i) <= 0 && scorebox(i) > tneg_score[index])) && sodbox(i) > t_sod[index] && scorebox(i) > scorebox_max ) {
                            scorebox_max = scorebox(i);
                            index_scorebox_max = index;
                            i_scorebox_max = i;
                            if(debug) std::cout << std::setprecision(3) << "step1 - flag=1? - index scorebox sodbox corep distance " << index << " " << scorebox(i) << " " << sodbox(i) << " " << corep(i) << " " << dist << std::endl;
                        } else {
                            update[i] = 1;
                            if(debug) std::cout << std::setprecision(3) << "step1 - flag=0? - index scorebox sodbox corep distance " << index << " " << scorebox(i) << " " << sodbox(i) << " " << corep(i) << " " << dist << std::endl;
                        }
                        count_updates++;
                    }
                }
            }
            
            // step 2. update the flags
            if( count_updates > 0) {
                // if all observations are good ones ...
                if( scorebox_max < 0) {
                    for(int i = 0; i < p_outer; i++) {
                        if ( update[i] == 1) {
                            int index = neighbour_indices[i];
                            score[index] = scorebox(i);
                            rep[index] = corep(i);
                            sod[index] = sodbox(i);
                            num_inner[index] = p_inner;
                            horizontal_scale[index] = Dhbox_mean;
                            an_inc[index] = ainc(i);
                            an_res[index] = ares(i);
                            cv_res[index] = cvres(i);
                            innov[index] = d(i);
                            idi[index] = idibox(i);
                            idiv[index] = idivbox(i);
                            sig2o[index] = sig2obox;
                            flags[index] = 0;
                            if(debug) {
                                int j = index;
                                std::cout << std::setprecision(3) << "    step2 - flag=0 - index innov ares ainc cvres idi idiv: " << j << " " << innov[j] << " " << an_res[j] << " " << an_inc[j] << " " << cv_res[j] << " " << idi[j] << " " << idiv[j] << " "  << std::endl;
                            }
                        }
                    }
                // ... otherwise, flag just the one with the highest score
                } else {
                    score[index_scorebox_max] = scorebox_max;
                    rep[index_scorebox_max] = corep(i_scorebox_max);
                    sod[index_scorebox_max] = sodbox(i_scorebox_max);
                    num_inner[index_scorebox_max] = p_inner;
                    horizontal_scale[index_scorebox_max] = Dhbox_mean;
                    an_inc[index_scorebox_max] = ainc(i_scorebox_max);
                    an_res[index_scorebox_max] = ares(i_scorebox_max);
                    cv_res[index_scorebox_max] = cvres(i_scorebox_max);
                    innov[index_scorebox_max] = d(i_scorebox_max);
                    idi[index_scorebox_max] = idibox(i_scorebox_max);
                    idiv[index_scorebox_max] = idivbox(i_scorebox_max);
                    sig2o[index_scorebox_max] = sig2obox;
                    flags[index_scorebox_max] = 1;
                    thrown_out++;
                    if(debug) {
                        int j = index_scorebox_max;
                        std::cout << std::setprecision(3) << "    step2 - flag=1 - index innov ares ainc cvres idi idiv: " << j << " " << innov[j] << " " << an_res[j] << " " << an_inc[j] << " " << cv_res[j] << " " << idi[j] << " " << idiv[j] << " "  << std::endl;
                    }
                }
            }
            count_oi++; // one more OI
        }  // end loop over observations
        if(thrown_out == 0) {
            if(iteration + 1 < num_iterations)
                std::cout << "Stopping early after " << iteration + 1<< " iterations" << std::endl;
            break;
        }
        std::cout << "Removing " << thrown_out << " Number of OI " << count_oi << std::endl;
        double e_time0 = titanlib::util::clock();
        std::cout << e_time0 - s_time0 << std::endl;
    } // end of SCT iterations

    return flags;
}
// end SCT //


vec compute_vertical_profile(const vec& elevs, const vec& oelevs, const vec& values, int num_min_prof, double min_elev_diff) {
    // Starting value guesses
    double gamma = -0.0065;
    double a = 5.0;

    double meanT = std::accumulate(values.begin(), values.end(), 0.0) / values.size();
    double exact_p10 = titanlib::util::compute_quantile(0.10, elevs);
    double exact_p90 = titanlib::util::compute_quantile(0.90, elevs);

    // optimize inputs for VP (using Nelder-Mead Simplex algorithm)
    const gsl_multimin_fminimizer_type *T = gsl_multimin_fminimizer_nmsimplex2;
    gsl_multimin_fminimizer *s = NULL;
    gsl_vector *ss;
    gsl_multimin_function vp_optim;

    int iter = 0;
    int status;
    double size;

    // data (params) that needs to be passed into vp
    double nd = (double) elevs.size(); // cast + resize to double
    double nod = (double) oelevs.size(); // cast + resize to double
    std::vector<double> delevs(elevs.begin(), elevs.end());
    std::vector<double> doelevs(oelevs.begin(), oelevs.end());
    std::vector<double> dvalues(values.begin(), values.end());
    double * dpelevs = delevs.data();
    double * dpoelevs = doelevs.data();
    double * dpvalues = dvalues.data();
    double * data[3] = {&nd, dpelevs, dpvalues};

    // Check if terrain is too flat
    double z05 = titanlib::util::compute_quantile(0.05, elevs);
    double z95 = titanlib::util::compute_quantile(0.95, elevs);

    // should we use the basic or more complicated vertical profile?
    bool use_basic = elevs.size() < num_min_prof || (z95 - z05) < min_elev_diff;

    gsl_vector* input;
    if(use_basic) {
        vp_optim.n = 2;
        input = gsl_vector_alloc(vp_optim.n);
        gsl_vector_set(input, 0, meanT);
        gsl_vector_set(input, 1, gamma);
        vp_optim.f = basic_vertical_profile_optimizer_function;
    }
    else {
        vp_optim.n = 5;
        input = gsl_vector_alloc(vp_optim.n);
        gsl_vector_set(input, 0, meanT);
        gsl_vector_set(input, 1, gamma);
        gsl_vector_set(input, 2, a);
        gsl_vector_set(input, 3, exact_p10);
        gsl_vector_set(input, 4, exact_p90);
        vp_optim.f = vertical_profile_optimizer_function;
    }
    ss = gsl_vector_alloc (vp_optim.n);

    gsl_vector_set_all (ss, 1.0);
    vp_optim.params = data;

    s = gsl_multimin_fminimizer_alloc (T, vp_optim.n);
    gsl_multimin_fminimizer_set (s, &vp_optim, input, ss);
    do {
        iter++;
        status = gsl_multimin_fminimizer_iterate(s);

        if (status)
            break;

        size = gsl_multimin_fminimizer_size (s);
        status = gsl_multimin_test_size (size, 1e-2);
    }
    while (status == GSL_CONTINUE && iter < 100);

    // then actually calculate the vertical profile using the minima
    vec vp;
    if(use_basic) {
        vp = basic_vertical_profile(nod, dpoelevs, gsl_vector_get(s->x, 0), gsl_vector_get(s->x, 1));
        // std::cout << "meanT=" << gsl_vector_get(s->x, 0) << " gamma=" << gsl_vector_get(s->x, 1) << std::endl;
    }
    else {
        // then actually calculate the vertical profile using the minima
        vp = vertical_profile(nod, dpoelevs,  gsl_vector_get(s->x, 0), gsl_vector_get(s->x, 1),
                gsl_vector_get(s->x, 2), gsl_vector_get(s->x, 3), gsl_vector_get(s->x, 4));
    }

    gsl_vector_free(input);
    gsl_vector_free(ss);
    gsl_multimin_fminimizer_free (s);

    return vp;
}

double basic_vertical_profile_optimizer_function(const gsl_vector *v, void *data) {
    double **p = (double **)data;
    int n = (int) *p[0]; // is of type double but should be an int
    double *dpelevs = p[1];
    double *dpvalues = p[2];

    // the parameters to mess with
    double meanT = gsl_vector_get(v,0);
    double gamma = gsl_vector_get(v,1);

    // give everything to vp to compute t_out
    vec t_out = basic_vertical_profile(n, dpelevs, meanT, gamma);
    // RMS
    float total = 0;
    for(int i=0; i<n; i++) {
        total += pow((t_out[i]-dpvalues[i]),2);
    }
    double value = log(pow((total / n),0.5));
    return value;
}

vec basic_vertical_profile(const int n, const double *elevs, const double t0, const double gamma) {
    vec t_out(n, -999);
    for(int i=0; i<n; i++)
        t_out[i] = t0 + gamma*elevs[i];
    return t_out;
}

double vertical_profile_optimizer_function(const gsl_vector *v, void *data) {
    double **p = (double **)data;
    int n = (int) *p[0]; // is of type double but should be an int
    double *dpelevs = p[1];
    double *dpvalues = p[2];

    // the parameters to mess with
    double meanT = gsl_vector_get(v,0);
    double gamma = gsl_vector_get(v,1);
    double a = gsl_vector_get(v,2);
    double exact_p10 = gsl_vector_get(v,3);
    double exact_p90 = gsl_vector_get(v,4);

    // give everything to vp to compute t_out
    vec t_out = vertical_profile(n, dpelevs, meanT, gamma, a, exact_p10, exact_p90);
    // RMS
    double total = 0;
    for(int i=0; i<n; i++) {
        total += pow((t_out[i]-dpvalues[i]),2);
    }
    double value = log(pow((total / n),0.5));
    return value;
}

vec vertical_profile(const int n, const double *elevs, const double t0, const double gamma, const double a, const double h0, const double h1i) {
    double h1 = h0 + fabs(h1i); // h1<-h0+abs(h1i)
    // loop over the array of elevations
    vec t_out;
    t_out.resize(n, -999);

    for(int i=0; i<n; i++) {
        // define some bools
        bool z_le_h0 = elevs[i] <= h0; // z.le.h0<-which(z<=h0)
        bool z_ge_h1 = elevs[i] >= h1; // z.ge.h1<-which(z>=h1)
        bool z_in = (elevs[i]>h0 && elevs[i]<h1); // z.in<-which(z>h0 & z<h1)
        if(z_le_h0) {
            t_out[i] = t0-gamma*elevs[i]-a;
        }
        if(z_ge_h1) {
            t_out[i] = t0-gamma*elevs[i];
        }
        if(z_in) {
            t_out[i] = t0-gamma*elevs[i]-a/2*(1+cos(M_PI*(elevs[i]-h0)/(h1-h0)));
        }
    }
    return t_out;
}


//----------------------------------------------------------------------------//
// HELPER FUNCTIONS

bool invert_matrix(const boost::numeric::ublas::matrix<float>& input, boost::numeric::ublas::matrix<float>& inverse) {
    using namespace boost::numeric::ublas;
    typedef permutation_matrix<std::size_t> pmatrix;
    // create a working copy of the input
    matrix<float> A(input);
    // create a permutation matrix for the LU-factorization
    pmatrix pm(A.size1());
    // perform LU-factorization
    int res = lu_factorize(A,pm);
    if( res != 0 ) return false;
    // create identity matrix of "inverse"
    inverse.assign(boost::numeric::ublas::identity_matrix<float>(A.size1()));
    // backsubstitute to get the inverse
    lu_substitute(A, pm, inverse);
    return true;
}


ivec remove_flagged(ivec indices, ivec flags) {
    ivec removed;
    removed.reserve(indices.size());
    for(int i=0; i<indices.size(); i++) {
        if(flags[indices[i]] != 1 ) {
            removed.push_back(indices[i]);
        }
    }
    return removed;
}
