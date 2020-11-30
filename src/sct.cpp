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
//ivec remove_flagged(ivec indices, ivec flags);
ivec remove_flagged(ivec indices, ivec flags, vec dist, vec &dist_updated, int keep);

bool oi_adhoc( const vec& lats, const vec& lons, const vec& elevs, const vec& values, const boost::numeric::ublas::vector<float>& d, float min_horizontal_scale, float max_horizontal_scale, int kth_close, float vertical_scale, float sig2o_min, float sig2o_max, const vec& eps2, bool debug, const float& na, double& Dh_mean, double& sig2o, double& sig2o_mean_var, double& score_mean, double& score_var, double& score_mean_var, boost::numeric::ublas::vector<float>& ainc, boost::numeric::ublas::vector<float>& ares, boost::numeric::ublas::vector<float>& cvres, boost::numeric::ublas::vector<float>& idi, boost::numeric::ublas::vector<float>& idiv, boost::numeric::ublas::vector<float>& score_numerator, boost::numeric::ublas::vector<float>& rep_numerator);

vec compute_vertical_profile(const vec& elevs, const vec& oelevs, const vec& values, int num_min_prof, double min_elev_diff, bool debug);

vec compute_vertical_profile_Theil_Sen(const vec& elevs, const vec& oelevs, const vec& values, int num_min_prof, double min_elev_diff, bool debug);

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
 + Coefficient of representativeness (rep, p_outer-vector)
   rep = innovation * analysis residuals / sig2o, LUS10 Eq.(32) numerator
 + Spatial Outlier Detection (SOD) score (sod, p_outer-vector)
    used to distinguish between acceptable REs and GEs
     sod = (deviation from SCT-score areal average)^2 / spread 

 Basic Algorithm: 
    Loop over the s observations (observation index is "curr"). Select the subset of p_outer closest observations to the curr-th observation, inside a predefined outer circle, then perform SCT considering only this subset.
    Define an inner circle around the curr-th observation and only observations inside the inner circle can be flagged.
    Observations are candidate "bad" observations if the corresponding scores are a) larger than pre-set thresholds AND b) outliers compared to the statistics of score inside the outer circle
      NOTE: condition b) outlier if sod>threshold. Helps to avoid flagging REs as GEs, at least in those regions where enough observations are avaialable
    Among all candidate bad observations, only the one with the larges score is flagged as "bad"
    If there are not candidate bad observations, then flag all the observations inside the inner circle as "good" ones
 
 The basic algorithm is repeated several times, each time learing from the previous iteration and testing only those observations that are left without a decision. Note that the first iteration is used only to identify bad observations, as a measure of precaution.
 The SCT-loop should end when 0 observations are flagged as "bad" ones.

 Finally, a round that test only bad observations using good observations is performed. A bad observation can be "saved" if it passes the final round.

 Returned values:
  flags. -999 = not checked; 0 = passed (good); 1 = failed (bad); 11 = isolated (<2 inside inner); 12 = isolated (<num_min inside outer)

*/
    double s_time = titanlib::util::clock();
    
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
    if(background_elab_type != "vertical_profile" && background_elab_type != "vertical_profile_Theil_Sen" && background_elab_type != "mean_outer_circle" && background_elab_type != "external")
        throw std::invalid_argument("background_elab_type must be one of vertical_profile, vertical_profile_Theil_Sen, mean_outer_circle or external");
    if(background_elab_type == "vertical_profile" && num_min_prof<0)
        throw std::invalid_argument("num_min_prof must be >=0");
    if( (obs_to_check.size() != p && obs_to_check.size() != 1 && obs_to_check.size() !=0) ) 
        throw std::invalid_argument("'obs_to_check' has an invalid length");
    if(background_elab_type == "external" &&  background_values.size() != p)
        throw std::runtime_error("Background vector dimension mismatch");

    // initializations
    float na = -999.; // code for "Not Available". Any type of missing data
    ivec flags( p, na);
    vec perc_flag0( p, na);
    bool set_zeros = false;
 
    score.clear();
    score.resize( p, na);
    sod.clear();
    sod.resize( p, na);
    num_inner.clear();
    num_inner.resize( p, na);
    horizontal_scale.clear();
    horizontal_scale.resize( p, na);
    an_inc.clear();
    an_inc.resize( p, na);
    an_res.clear();
    an_res.resize( p, na);
    cv_res.clear();
    cv_res.resize( p, na);
    innov.clear();
    innov.resize( p, na);
    idi.clear();
    idi.resize( p, na);
    idiv.clear();
    idiv.resize( p, na);
    sig2o.clear();
    sig2o.resize( p, na);
    rep.clear();
    rep.resize( p, na);

    /* thresholds used to decide if the innovations are so small that the sct can be skipped
       since the observations stay very close to the background */
    vec innov_small( p, na);
    for(int j=0; j < p; j++) {
        innov_small[j] = std::pow( sig2o_min + sig2o_min / eps2[j], 0.5);
    }


    /* preliminary range check
       note: the sct can mistake obvious large errors for good observations 
             if they happen to be close enough that they support each other */
    int thrown_out = 0; 
    for(int j=0; j < p; j++) {
        if((values[j] < value_min) || (values[j] > value_max)) { 
          flags[j] = 1;
          thrown_out++;
        }
    }
    if(debug) std::cout << "=============================================== " << p << std::endl;
    if(debug) std::cout << "Number of observations to test is " << p << std::endl;
    if(debug) std::cout << "Range check is removing " << thrown_out << " observations " << std::endl;

    // if obs_to_check is empty then check all
    bool check_all = obs_to_check.size() != p;

    // KDtree has to do with fast computation of distances
    titanlib::KDTree tree(lats, lons);


    // SCT iterations
    for(int iteration = 0; iteration < num_iterations; iteration++) {

        if(debug) std::cout << " +++++ Iteration " << iteration << " ++++++++++++++++" << std::endl;

        double s_time0 = titanlib::util::clock();
        
        // reset this number each loop (this is for breaking if we don't throw anything new out)
        int thrown_out = 0; 
        
        // diagnostic. count the number of times OI is performed
        int count_oi = 0;

        // loop over observations
        for(int curr=0; curr < p; curr++) {

            if(debug) std::cout << "===> curr " << curr << " ===============" << std::endl;

            // observation not to be checked
            if(!check_all && obs_to_check[curr] != 1) {
                if(debug) std::cout << "..observation not to check " << curr << std::endl;
                continue;
            }

            // break out if station already flagged
            if( flags[curr] >= 0) {
                if(debug) std::cout << "..checked " << curr << std::endl;
                continue;
            }

            // if no neighbours inside the inner circle, then distinction between RE and GE not reliable
            vec distances_inner, distances_inner_tmp;
            ivec neighbour_indices_inner = tree.get_neighbours_with_distance(lats[curr], lons[curr], inner_radius, num_max, true, distances_inner_tmp);
            neighbour_indices_inner = remove_flagged( neighbour_indices_inner, flags, distances_inner_tmp, distances_inner, curr);
            if( neighbour_indices_inner.size() < 2) {
                flags[curr] = 11;
                if(debug) std::cout << "@@isolated (inner) " << curr << std::endl;
                continue;
            }
            int p_inner = neighbour_indices_inner.size();
 
            // get all neighbours that are close enough (inside outer circle)
            vec distances;
            vec distances_tmp;
            ivec neighbour_indices = tree.get_neighbours_with_distance(lats[curr], lons[curr], outer_radius, num_max, true, distances_tmp);
            // curr-th observation is in the last position of the vector
            neighbour_indices = remove_flagged( neighbour_indices, flags, distances_tmp, distances, curr);
            if(neighbour_indices.size() < num_min) {
                flags[curr] = 12;
                if(debug) std::cout << "@@isolated (outer) " << curr << std::endl;
                continue; 
            }

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
                // vertical profile is independent from the curr-th observation
                std::vector<float> elevs_box1(elevs_box.begin(), elevs_box.end() - 1);
                std::vector<float> values_box1(values_box.begin(), values_box.end() - 1);
                bvalues_box = compute_vertical_profile(elevs_box1, elevs_box, values_box1, num_min_prof, min_elev_diff, debug);
            } else if( background_elab_type == "vertical_profile_Theil_Sen") {
                // vertical profile is independent from the curr-th observation
                std::vector<float> elevs_box1(elevs_box.begin(), elevs_box.end() - 1);
                std::vector<float> values_box1(values_box.begin(), values_box.end() - 1);
                bvalues_box = compute_vertical_profile_Theil_Sen(elevs_box1, elevs_box, values_box1, num_min_prof, min_elev_diff, debug);
            } else if( background_elab_type == "mean_outer_circle"){
                std::vector<float> values_box1(values_box.begin(), values_box.end() - 1);
                double mean_val = std::accumulate(values_box1.begin(), values_box1.end(), 0.0) / ( p_outer - 1);
                bvalues_box.resize(p_outer, mean_val);
            } else if( background_elab_type == "external"){
                bvalues_box = titanlib::util::subset(background_values, neighbour_indices);
            } 
            if(debug) std::cout << "... background ok ..." << std::endl;
            
            // Compute innovations (= observations - background)
            boost::numeric::ublas::vector<float> innov_box(p_outer);
            bool small_innov = true;
            for(int i=0; i<p_outer; i++) {
                int index = neighbour_indices[i];
                // check the background is in a range of acceptable values
                if(bvalues_box[i] < value_min) bvalues_box[i] = value_min;
                if(bvalues_box[i] > value_max) bvalues_box[i] = value_max;
                // innovations
                innov_box(i) = values_box[i] - bvalues_box[i];
                if ( fabs( innov_box(i)) > innov_small[index]) small_innov = false;
                if(debug) std::cout << std::setprecision(3) << " backg - index elev obs backg " << neighbour_indices[i] << " " << elevs_box[i] << " " << values_box[i] << " " << bvalues_box[i] << std::endl;
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
                        if (dist <= inner_radius) {
                            score[index] = 0;
                            rep[index] = 0;
                            sod[index] = 0;
                            num_inner[index] = p_inner;
                            innov[index] = innov_box(i);
                            horizontal_scale[index] = na;
                            an_inc[index] = na;
                            an_res[index] = na;
                            cv_res[index] = na;
                            idi[index] = na;
                            idiv[index] = na;
                            sig2o[index] = na;
                            flags[index] = 0;
                            if (debug) std::cout << " small_innov - index " << index << std::endl;
                        }
                    }
                }
                continue; // jump to the next observation
            }

            // Optimal Interpolation, ad-hoc implementation for SCT
            double Dhbox_mean, sig2obox, sig2obox_mean_var, scorebox_mean, scorebox_var, scorebox_mean_var; 
            boost::numeric::ublas::vector<float> an_inc_box(p_outer), an_res_box(p_outer), cv_res_box(p_outer), rep_box_numerator(p_outer), idi_box(p_outer), idiv_box(p_outer), scorebox_numerator(p_outer);
            bool res = oi_adhoc( lats_box, lons_box, elevs_box, values_box, innov_box, min_horizontal_scale, max_horizontal_scale, kth_closest_obs_horizontal_scale, vertical_scale, sig2o_min, sig2o_max, eps2_box, debug, na, Dhbox_mean, sig2obox, sig2obox_mean_var, scorebox_mean, scorebox_var, scorebox_mean_var, an_inc_box, an_res_box, cv_res_box, idi_box, idiv_box, scorebox_numerator, rep_box_numerator);
            // this can only happen with problems during the matrix inversion, which should not happen...
            if ( !res) {
                flags[curr] = 100;
                if(debug) std::cout << " oi - flags=100 - index " << curr << std::endl; 
                continue;
            }
            if(debug) {
                for(int i=0; i < p_outer; i++) {
                    std::cout << std::setprecision(3) << " oi - index elev obs backg an cv_an idi idiv: " << neighbour_indices[i] << " " << elevs_box[i] << " " << values_box[i] << " " << values_box[i]-innov_box(i) << " " << values_box[i]-an_res_box(i) << " " << values_box[i]-cv_res_box(i) << " " << idi_box(i) << " " << idiv_box(i) << " "  << std::endl;
                }
            }
            
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
            int count_flag0 = 0;
            ivec update(p_outer, 0);  // Keep track of which observations need an update (for flags=0)
            boost::numeric::ublas::vector<float> scorebox(p_outer), sodbox(p_outer), rep_box(p_outer);
            for(int i = 0; i < p_outer; i++) {
                int index = neighbour_indices[i];
                if ( flags[index] == 0) count_flag0++; // number of good observations in the outer circle
                if ( flags[index] == na) {
                    float dist = distances[i];
                    if(dist <= inner_radius) {
//                        scorebox(i) = scorebox_numerator(i) / ( sig2obox + std::sqrt(sig2obox_mean_var));
                        scorebox(i) = scorebox_numerator(i) / sig2obox;
                        sodbox(i) = std::pow( ( scorebox(i) - scorebox_mean), 2) / (scorebox_var + scorebox_mean_var);
                        rep_box(i) = rep_box_numerator(i) / sig2obox;
                        // condition defining bad observations 
                        if( ((cv_res_box(i) > 0 && scorebox(i) > tpos_score[index]) || (cv_res_box(i) <= 0 && scorebox(i) > tneg_score[index])) && sodbox(i) > t_sod[index] && scorebox(i) > scorebox_max ) {
                            // candidate for bad observation
                            scorebox_max = scorebox(i);
                            index_scorebox_max = index;
                            i_scorebox_max = i;
                            if(debug) std::cout << std::setprecision(3) << " decision - step1 - flag=1? - index scorebox sodbox corep distance " << index << " " << scorebox(i) << " " << sodbox(i) << " " << rep_box(i) << " " << dist << std::endl;
                        } else {
                            // candidate for good observation
                            update[i] = 1;
                            if(debug) std::cout << std::setprecision(3) << " decision - step1 - flag=0? - index scorebox sodbox corep distance " << index << " " << scorebox(i) << " " << sodbox(i) << " " << rep_box(i) << " " << dist << std::endl;
//                            if ( scorebox(i) > score[index]) {
//                                // update diagnostics, note that we keep the "worst" scores
//                                score[index] = scorebox(i);
//                                rep[index] = rep_box(i);
//                                sod[index] = sodbox(i);
//                                num_inner[index] = p_inner;
//                                horizontal_scale[index] = Dhbox_mean;
//                                an_inc[index] = an_inc_box(i);
//                                an_res[index] = an_res_box(i);
//                                cv_res[index] = cv_res_box(i);
//                                innov[index] = innov_box(i);
//                                idi[index] = idi_box(i);
//                                idiv[index] = idiv_box(i);
//                                sig2o[index] = sig2obox;
//                            }
                        }
                        count_updates++;
                    }
                }
            }
            /* percentage of good observations used to decide the flags
               the more good observations we use, the more we are confident in flagging a bad one 
               note: this parameter is used to optimize the running time in the final decision on 
               bad observations */
            float perc_flag0_curr = count_flag0 / p_outer;
           
            // step 2. update the flags
            if( count_updates > 0) {
                /* if all observations in the inner circle are good ones, then flag them with
                   flag=0. Do it from the 2nd iteration onward. The first iteration is only
                   used to remove the worst data from the dataset
                */
                if( scorebox_max < 0) {
                    for(int i = 0; i < p_outer; i++) {
                        if ( update[i] == 1) {
                            int index = neighbour_indices[i];
                            score[index] = scorebox(i);
                            rep[index] = rep_box(i);
                            sod[index] = sodbox(i);
                            num_inner[index] = p_inner;
                            horizontal_scale[index] = Dhbox_mean;
                            an_inc[index] = an_inc_box(i);
                            an_res[index] = an_res_box(i);
                            cv_res[index] = cv_res_box(i);
                            innov[index] = innov_box(i);
                            idi[index] = idi_box(i);
                            idiv[index] = idiv_box(i);
                            sig2o[index] = sig2obox;
                            if ( iteration > 0) flags[index] = 0;
                            if(debug) 
                                std::cout << " decision - step2 - flag=0 - index " << index << " " << std::endl;
                        }
                    }
                /* ... otherwise, flag as bad just the one with the highest score
                       in fact, the sct assumes all the observations to be good ones
                       then a bad observations affect the statistics of all the others.
                       note: the first iteration can only set flag=1 (not flag=0) */
                } else {
                    score[index_scorebox_max] = scorebox_max;
                    rep[index_scorebox_max] = rep_box(i_scorebox_max);
                    sod[index_scorebox_max] = sodbox(i_scorebox_max);
                    num_inner[index_scorebox_max] = p_inner;
                    horizontal_scale[index_scorebox_max] = Dhbox_mean;
                    an_inc[index_scorebox_max] = an_inc_box(i_scorebox_max);
                    an_res[index_scorebox_max] = an_res_box(i_scorebox_max);
                    cv_res[index_scorebox_max] = cv_res_box(i_scorebox_max);
                    innov[index_scorebox_max] = innov_box(i_scorebox_max);
                    idi[index_scorebox_max] = idi_box(i_scorebox_max);
                    idiv[index_scorebox_max] = idiv_box(i_scorebox_max);
                    sig2o[index_scorebox_max] = sig2obox;
                    perc_flag0[index_scorebox_max] = perc_flag0_curr;
                    flags[index_scorebox_max] = 1;
                    thrown_out++;
                    if(debug) 
                        std::cout << " decision - step2 - flag=1 - index " << index_scorebox_max << " " << std::endl;
                }
            }
            count_oi++; // one more OI
        }  // end loop over observations
        std::cout << "Removing " << thrown_out << " observations. Number of OI " << count_oi << std::endl;
        double e_time0 = titanlib::util::clock();
        std::cout << e_time0 - s_time0 << " secs" << std::endl;
        if(thrown_out == 0) {
            if ( iteration == 0) set_zeros = true;
            if(iteration + 1 < num_iterations)
                std::cout << "Stopping early after " << iteration + 1<< " iterations" << std::endl;
            break;
        }
    } // end of SCT iterations

    /* enter if the first iteration has been completed without flagging any bad observation
       the loop flags all the observations as good ones. This is prevented in the SCT 
       iterations. */
    if ( set_zeros) {
        int count_zero = 0;
        for(int curr=0; curr < p; curr++) {
          if ( flags[curr] == na && score[curr] > na) {
            flags[curr] = 0;
            count_zero++;
          }
        }
        if(debug) std::cout << "Remaining " << count_zero << " observations set to \"good\" " << std::endl;
    }

    /* final check on the bad observations
       it may happen that a good observation is flagged as a bad one because of the order
       the SCT has been done and the uncertainty made in the estimation of the background 
       (remember that bad observations may have been used to get the background).
       This final check is made only on bad observations and uses only good observations. */
    if(debug) std::cout << " +++++ Final check on bad observations ++++++++++++++++" << std::endl;
    double s_time0 = titanlib::util::clock();
    int rethrown_out = 0; 
    int saved = 0; 
    int count_oi = 0;
    // loop over observations
    for(int curr=0; curr < p; curr++) {

        if(debug) std::cout << "---> curr " << curr << "---------------" << std::endl;
        
        if((values[curr] < value_min) || (values[curr] > value_max)) {
            if(debug) std::cout << "..range check " << curr << std::endl;
            rethrown_out++; 
            continue;
        }

        // observation not to be checked
        if(!check_all && obs_to_check[curr] != 1) {
            if(debug) std::cout << "..not to check " << curr << std::endl;
            continue;
        }

        /* break out if station not flagged as bad or if flagged as bad
           when more than 50% of the observations used was flagged as good */
        if( flags[curr] != 1 || ( flags[curr] == 1 && perc_flag0[curr] > 0.5)) {
            if(debug) std::cout << "... flag!=1 or we trust the previous decision flag=1 " << curr << std::endl;
            continue;
        }

        // if no neighbours inside the inner circle, then distinction between RE and GE not reliable
        vec distances_inner;
        vec distances_inner_tmp;
        ivec neighbour_indices_inner = tree.get_neighbours_with_distance(lats[curr], lons[curr], inner_radius, num_max, true, distances_inner_tmp);
        neighbour_indices_inner = remove_flagged(neighbour_indices_inner, flags, distances_inner_tmp, distances_inner, curr);
        if(neighbour_indices_inner.size() < 2) {
            flags[curr] = 11;
            if(debug) std::cout << "@@isolated (inner) " << curr << std::endl;
            continue;
        }
        int p_inner = neighbour_indices_inner.size();

        // get all neighbours that are close enough (inside outer circle)
        vec distances;
        vec distances_tmp;
        ivec neighbour_indices = tree.get_neighbours_with_distance(lats[curr], lons[curr], outer_radius, num_max, true, distances_tmp);
        // curr-th observation is in the last position of the vector
        neighbour_indices = remove_flagged(neighbour_indices, flags, distances_tmp, distances, curr);
        if(neighbour_indices.size() < num_min) {
            flags[curr] = 12;
            if(debug) std::cout << "@@isolated (outer) " << curr << std::endl;
            continue; 
        }

        // call SCT with this box(=outer circle) 
        vec lons_box = titanlib::util::subset(lons, neighbour_indices);
        vec elevs_box = titanlib::util::subset(elevs, neighbour_indices);
        vec lats_box = titanlib::util::subset(lats, neighbour_indices);
        vec values_box = titanlib::util::subset(values, neighbour_indices);
        vec eps2_box = titanlib::util::subset(eps2, neighbour_indices);
        int p_outer = neighbour_indices.size();
        int i_curr = p_outer - 1; // pointer to curr-th observations
        if(debug) std::cout << "p_outer p_inner " << p_outer << " " << p_inner << " " << std::endl;
        if(debug) {
            int index = neighbour_indices[i_curr];
            std::cout << "i_curr index " << i_curr << " " << index << std::endl;
        }

        // Compute the background, if needed
        vec bvalues_box;
        if( background_elab_type == "vertical_profile") {
            // vertical profile is independent from the curr-th observation
            std::vector<float> elevs_box1(elevs_box.begin(), elevs_box.end() - 1);
            std::vector<float> values_box1(values_box.begin(), values_box.end() - 1);
            bvalues_box = compute_vertical_profile(elevs_box1, elevs_box, values_box1, num_min_prof, min_elev_diff, debug);
        } else if( background_elab_type == "vertical_profile_Theil_Sen") {
            // vertical profile is independent from the curr-th observation
            std::vector<float> elevs_box1(elevs_box.begin(), elevs_box.end() - 1);
            std::vector<float> values_box1(values_box.begin(), values_box.end() - 1);
            bvalues_box = compute_vertical_profile_Theil_Sen(elevs_box1, elevs_box, values_box1, num_min_prof, min_elev_diff, debug);
        } else if( background_elab_type == "mean_outer_circle"){
            std::vector<float> values_box1(values_box.begin(), values_box.end() - 1);
            double mean_val = std::accumulate(values_box.begin(), values_box.end(), 0.0) / values_box.size();
            bvalues_box.resize(p_outer, mean_val);
        } else if( background_elab_type == "external"){
            bvalues_box = titanlib::util::subset(background_values, neighbour_indices);
        } 
        if(debug) std::cout << "... background ok ..." << std::endl;
        
        // Compute innovations (= observations - background)
        boost::numeric::ublas::vector<float> innov_box(p_outer);
        bool small_innov = true;
        for(int i=0; i<p_outer; i++) {
            int index = neighbour_indices[i];
            // check the background is in a range of acceptable values
            if(bvalues_box[i] < value_min) bvalues_box[i] = value_min;
            if(bvalues_box[i] > value_max) bvalues_box[i] = value_max;
            // innovations
            innov_box(i) = values_box[i] - bvalues_box[i];
            if ( fabs( innov_box(i)) > innov_small[index]) small_innov = false;
            if(debug) std::cout << std::setprecision(3) << " backg - index elev obs backg " << neighbour_indices[i] << " " << elevs_box[i] << " " << values_box[i] << " " << bvalues_box[i] << std::endl;
        }
        
        /* if observations and background are almost identical, then take a shortcut
           flag = 0 for all observations in the inner circle
           NOTE: when innovations are all exectly =0, then the routine crashes
         */
        if (small_innov) {
            score[curr] = 0;
            rep[curr] = 0;
            sod[curr] = 0;
            num_inner[curr] = p_inner;
            innov[curr] = innov_box(i_curr);
            horizontal_scale[curr] = na;
            an_inc[curr] = na;
            an_res[curr] = na;
            cv_res[curr] = na;
            idi[curr] = na;
            idiv[curr] = na;
            sig2o[curr] = na;
            flags[curr] = 0;
            saved++;
            if (debug) std::cout << " saved - small_innov - index " << curr << std::endl;
            continue; // jump to the next observation
        }
        // Optimal Interpolation, ad-hoc implementation for SCT
        double Dhbox_mean, sig2obox, sig2obox_mean_var, scorebox_mean, scorebox_var, scorebox_mean_var; 
        boost::numeric::ublas::vector<float> an_inc_box(p_outer), an_res_box(p_outer), cv_res_box(p_outer), rep_box_numerator(p_outer), idi_box(p_outer), idiv_box(p_outer), scorebox_numerator(p_outer);
        bool res = oi_adhoc( lats_box, lons_box, elevs_box, values_box, innov_box, min_horizontal_scale, max_horizontal_scale, kth_closest_obs_horizontal_scale, vertical_scale, sig2o_min, sig2o_max, eps2_box, debug, na, Dhbox_mean, sig2obox, sig2obox_mean_var, scorebox_mean, scorebox_var, scorebox_mean_var, an_inc_box, an_res_box, cv_res_box, idi_box, idiv_box, scorebox_numerator, rep_box_numerator);
        // this can only happen with problems during the matrix inversion, which should not happen...
        if ( !res) {
            flags[curr] = 100;
            if(debug) std::cout << " oi - flags=100 - index " << curr << std::endl; 
            continue;
        }
        if(debug) {
            for(int i=0; i < p_outer; i++) {
                std::cout << std::setprecision(3) << " oi - index elev obs backg an cv_an idi idiv: " << neighbour_indices[i] << " " << elevs_box[i] << " " << values_box[i] << " " << values_box[i]-innov_box(i) << " " << values_box[i]-an_res_box(i) << " " << values_box[i]-cv_res_box(i) << " " << idi_box(i) << " " << idiv_box(i) << " "  << std::endl;
            }
        }

 
        // update flags inside the inner circle
//        double scorebox = scorebox_numerator(i_curr) / ( sig2obox + std::sqrt(sig2obox_mean_var));
        double scorebox = scorebox_numerator(i_curr) / sig2obox;
        double sodbox = std::pow( ( scorebox - scorebox_mean), 2) / (scorebox_var + scorebox_mean_var);
        double rep_box = rep_box_numerator(i_curr) / sig2obox;
        score[curr] = scorebox;
        rep[curr] = rep_box;
        sod[curr] = sodbox;
        num_inner[curr] = p_inner;
        horizontal_scale[curr] = Dhbox_mean;
        an_inc[curr] = an_inc_box(i_curr);
        an_res[curr] = an_res_box(i_curr);
        cv_res[curr] = cv_res_box(i_curr);
        innov[curr] = innov_box(i_curr);
        idi[curr] = idi_box(i_curr);
        idiv[curr] = idiv_box(i_curr);
        sig2o[curr] = sig2obox;
        // condition defining bad observations 
        if( ( ( cv_res_box(i_curr) > 0 && scorebox > tpos_score[curr]) || ( cv_res_box(i_curr) <= 0 && scorebox > tneg_score[curr])) && sodbox > t_sod[curr] ) {
            flags[curr] = 1;
            rethrown_out++;
            if(debug) std::cout << std::setprecision(3) << " rethrown out - flag=1 - index scorebox sodbox corep distance " << curr << " " << scorebox << " " << sodbox << " " << rep_box << " " << distances[i_curr] << std::endl;
        } else {
            flags[curr] = 0;
            saved++;
            if(debug) std::cout << std::setprecision(3) << " saved - flag=0 - index scorebox sodbox corep distance " << curr << " " << scorebox << " " << sodbox << " " << rep_box << " " << distances[i_curr] << std::endl;
        }
        count_oi++; // one more OI
    }  // end loop over observations


    std::cout << "Final decision. Thrown out " << rethrown_out << " Saved " << saved << " Number of OI " << count_oi << std::endl;
    double e_time0 = titanlib::util::clock();
    std::cout << e_time0 - s_time0 << "secs" << std::endl;
    
    //
    if(debug) {
        // loop over observations
        for(int curr=0; curr < p; curr++) {
            if( flags[curr] == 1) {
                double aux = an_res[curr] * cv_res[curr] / sig2o[curr];
//                if ( score[curr] < tpos_score[curr]) {
                if ( aux < tpos_score[curr]) {
                    std::cout << std::setprecision(3) << "@@BAD - curr an_res cv_res sig2o score aux: " << curr << " " << an_res[curr] << " " << cv_res[curr] << " " << sig2o[curr] << " " << score[curr] << " " << aux << std::endl; 
                } else {
                    std::cout << std::setprecision(3) << "  BAD - curr an_res cv_res sig2o score aux: " << curr << " " << an_res[curr] << " " << cv_res[curr] << " " << sig2o[curr] << " " << score[curr] << " " << aux << std::endl;
                }
            } else {
                double aux = an_res[curr] * cv_res[curr] / sig2o[curr];
//                if ( score[curr] < tpos_score[curr]) {
                if ( aux < tpos_score[curr]) {
                    std::cout << std::setprecision(3) << "   GOOD - curr an_res cv_res sig2o score aux: " << curr << " " << an_res[curr] << " " << cv_res[curr] << " " << sig2o[curr] << " " << score[curr] << " " << aux << std::endl;
                } else {
                    std::cout << std::setprecision(3) << " @@GOOD - curr an_res cv_res sig2o score aux: " << curr << " " << an_res[curr] << " " << cv_res[curr] << " " << sig2o[curr] << " " << score[curr] << " " << aux << std::endl;
                }
            }
        }
    }
    //
    
    std::cout << ">> Total Time " << e_time0 - s_time << "secs" << std::endl;

    //
    return flags;
}
// end SCT //

vec compute_vertical_profile_Theil_Sen(const vec& elevs, const vec& oelevs, const vec& values, int num_min_prof, double min_elev_diff, bool debug) {
    // Starting value guesses
    double gamma = -0.0065;

    double meanT = std::accumulate(values.begin(), values.end(), 0.0) / values.size();

    // Special case when all observations have the same elevation
    if ( std::min_element(elevs.begin(),elevs.end()) == std::max_element(elevs.begin(),elevs.end())) {
        vec vp( oelevs.size(), meanT);
        return vp;
    }

    // Check if terrain is too flat
    double z05 = titanlib::util::compute_quantile(0.05, elevs);
    double z95 = titanlib::util::compute_quantile(0.95, elevs);

    // should we use the basic or more complicated vertical profile?
    bool use_basic = elevs.size() < num_min_prof || (z95 - z05) < min_elev_diff;

    // Theil-Sen (Median-slope) Regression (Wilks (2019), p. 284)
    int n = elevs.size();
    int no = oelevs.size();
    int nm = n * (n - 1) / 2;
    std::vector<float> q(n), m(nm);
    double m_median;
    if ( use_basic) {
        m_median = gamma;
    } else {
        int k = 0;
        for(int i=0; i<(n-1); i++) {
            for(int j=(i+1); j<n; j++) {
                float dz = fabs( elevs[i] - elevs[j]);
                if ( dz < 1 ) {
                    m[k] = 0;
                } else {
                    m[k] = ( values[i] - values[j]) / ( elevs[i] - elevs[j]);
                }
                k++;
            }
        }
        m_median = titanlib::util::compute_quantile( 0.5, m);
    }
    for(int i=0; i<n; i++) {
        q[i] = values[i] - m_median * elevs[i];
    }
    double q_median = titanlib::util::compute_quantile( 0.5, q);
    vec vp(no);
    if (debug) std::cout << "Theil-Sen vp - m q " << m_median  << " " << q_median << std::endl;
    for(int i=0; i<no; i++) {
      vp[i] = q_median + m_median * oelevs[i];
    }
    //
    return vp;
}

vec compute_vertical_profile(const vec& elevs, const vec& oelevs, const vec& values, int num_min_prof, double min_elev_diff, bool debug) {
    // Starting value guesses
    double gamma = -0.0065;
    double a = 5.0;

    double meanT = std::accumulate(values.begin(), values.end(), 0.0) / values.size();

    // Special case when all observations have the same elevation
    if ( std::min_element(elevs.begin(),elevs.end()) == std::max_element(elevs.begin(),elevs.end())) {
        vec vp( oelevs.size(), meanT);
        return vp;
    }
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
        if(debug) std::cout << "vp - meanT=" << gsl_vector_get(s->x, 0) << " gamma=" << gsl_vector_get(s->x, 1) << std::endl;
    }
    else {
        // then actually calculate the vertical profile using the minima
        vp = vertical_profile(nod, dpoelevs,  gsl_vector_get(s->x, 0), gsl_vector_get(s->x, 1),
                gsl_vector_get(s->x, 2), gsl_vector_get(s->x, 3), gsl_vector_get(s->x, 4));
        if(debug) std::cout << "vp - t0 " << gsl_vector_get(s->x, 0) << " gamma=" << gsl_vector_get(s->x, 1) << " a=" << gsl_vector_get(s->x, 2) << " h0=" << gsl_vector_get(s->x, 3) << " h1i=" << gsl_vector_get(s->x, 4) << std::endl;
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
            t_out[i] = t0+gamma*elevs[i]-a;
        }
        if(z_ge_h1) {
            t_out[i] = t0+gamma*elevs[i];
        }
        if(z_in) {
            t_out[i] = t0+gamma*elevs[i]-a/2*(1+cos(M_PI*(elevs[i]-h0)/(h1-h0)));
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


ivec remove_flagged(ivec indices, ivec flags, vec dist, vec &dist_updated, int keep) {
    ivec removed;
    removed.reserve(indices.size());
    int j = -1;
    for(int i=0; i<indices.size(); i++) {
        if(flags[indices[i]] != 1 && indices[i] != keep) {
            removed.push_back(indices[i]);
            dist_updated.push_back(dist[i]);
        } else if ( indices[i] == keep) {
           j = i;
        }
    }
    if ( j >= 0) {
        removed.push_back(indices[j]);
        dist_updated.push_back(dist[j]);
    }
    return removed;
}

bool oi_adhoc( const vec& lats, const vec& lons, const vec& elevs, const vec& values, const boost::numeric::ublas::vector<float>& d, float min_horizontal_scale, float max_horizontal_scale, int kth_close, float vertical_scale, float sig2o_min, float sig2o_max, const vec& eps2, bool debug, const float& na, double& Dh_mean, double& sig2o, double& sig2o_mean_var, double& score_mean, double& score_var, double& score_mean_var, boost::numeric::ublas::vector<float>& ainc, boost::numeric::ublas::vector<float>& ares, boost::numeric::ublas::vector<float>& cvres, boost::numeric::ublas::vector<float>& idi, boost::numeric::ublas::vector<float>& idiv, boost::numeric::ublas::vector<float>& score_numerator, boost::numeric::ublas::vector<float>& rep_numerator) {
    using namespace boost::numeric::ublas;
    // init
    const int p = values.size();
    ainc.clear();
    ainc.resize( p, na);
    ares.clear();
    ares.resize( p, na);
    cvres.clear();
    cvres.resize( p, na);
    idi.clear();
    idi.resize( p, na);
    idiv.clear();
    idiv.resize( p, na);
    score_numerator.clear();
    score_numerator.resize( p, na);
    rep_numerator.clear();
    rep_numerator.resize( p, na);

    /* Compute Dh. 
     The location-dependent horizontal de-correlation lenght scale 
     used for the background error correlation matrix */
    boost::numeric::ublas::matrix<float> disth(p, p);
    boost::numeric::ublas::matrix<float> distz(p, p);
    boost::numeric::ublas::vector<float> Dh(p);

    for(int i=0; i < p; i++) {
        vec Dh_vector(p);
        for(int j=0; j < p; j++) {
            disth(i, j) = titanlib::util::calc_distance( lats[i], lons[i], lats[j], lons[j]);
            distz(i, j) = fabs( elevs[i] - elevs[j]);
            if(i != j) {
                if(i < j)
                    Dh_vector[j - 1] = disth(i, j);
                else if(i > j)
                    Dh_vector[j] = disth(i, j);
            }
        }
        Dh(i) = titanlib::util::findKclosest( kth_close, Dh_vector); // k-th closest observations
    }

    Dh_mean = std::accumulate(std::begin(Dh), std::end(Dh), 0.0) / Dh.size();
    if(Dh_mean < min_horizontal_scale) {
        if(debug) std::cout << "Dh_mean (<min_horizontal_scale) " << Dh_mean << std::endl;
        Dh_mean = min_horizontal_scale;
    }
    if(Dh_mean > max_horizontal_scale) {
        if(debug) std::cout << "Dh_mean (>max_horizontal_scale) " << Dh_mean << std::endl;
        Dh_mean = max_horizontal_scale;
    }
    if(debug) std::cout << "Dh_mean " << Dh_mean << std::endl;
    
    // Compute S + eps2*I and store it in S 
    boost::numeric::ublas::matrix<float> S(p,p);
    boost::numeric::ublas::matrix<float> Sinv(p,p);
    for(int i=0; i < p; i++) {
        for(int j=0; j < p; j++) {
            double value = std::exp(-.5 * std::pow((disth(i, j) / Dh_mean), 2) - .5 * std::pow((distz(i, j) / vertical_scale), 2));
            if(i==j) { // weight the diagonal, this also ensure an invertible matrix
                value = value + eps2[i];
            }
            S(i,j) = value;
        }
    }
    
    // Compute ( S + eps2 * I )^(-1) ans store it in Sinv
    bool b = invert_matrix(S, Sinv);
    if( !b)  return(false);

    // Definitions
    boost::numeric::ublas::vector<float> Zinv(p), Sinv_d(p), Sinv_1(p);

    // Matrix multiplications
    for(int i=0; i<p; i++) {
        S(i,i) -= eps2[i];
        double acc = 0;
        double acc1 = 0;
        for(int j=0; j<p; j++) {
            acc += Sinv(i,j)*d(j);
            acc1 += Sinv(i,j);
        }
        Sinv_d(i) = acc;
        Sinv_1(i) = acc1;
    }
   
    double M2a = 0;
    double M2b = 0;
    score_mean = 0;      // score mean
    sig2o = 0;           // sig2o mean
    for(int i=0; i<p; i++) {
        double acc = 0;
        double acc1 = 0;
        for(int j=0; j<p; j++) {
            acc += S(i,j)*Sinv_d(j);
            acc1 += S(i,j)*Sinv_1(j);
        }
        idi(i) = acc1;
        ainc(i) = acc;
        Zinv(i) = 1 / Sinv(i,i); 
        ares(i) = d(i) - ainc(i);
        cvres(i) = Zinv(i) * Sinv_d(i); 
        idiv(i) = 1 - Zinv(i) * Sinv_1(i); 
        // (Welford's online algorithm for mean and variance)
        // estimation of average score inside outer circle
        score_numerator(i) = cvres(i) * ares(i);
        float delta = score_numerator(i) - score_mean;
        score_mean += delta / (i+1);
        float delta2 = score_numerator(i) - score_mean;
        M2a += delta * delta2;
        // estimation of sig2o (observation error variance) 
        rep_numerator(i) = fabs( d(i) * ares(i)); // I've added fabs(...) 
        delta = rep_numerator(i) - sig2o;
        sig2o += delta / (i+1);
        delta2 = rep_numerator(i) - sig2o;
        M2b += delta * delta2;
    }

    // finalize estimation of sig2o sample variance
    double sig2o_var  = M2b / (p-1); 
    // sig2o variance of mean (see e.g. Taylor "an Intro to error analysis..." 1982, p. 102)
    sig2o_mean_var = sig2o_var / p;
    // one should have an idea of the range of allowed values for sig2o
    if(sig2o < sig2o_min) sig2o = sig2o_min;
    if(sig2o > sig2o_max) {
      sig2o = sig2o_max;
      sig2o_mean_var = 0; // high enough uncertainty is included in sig2o_max
    }
    if(debug) std::cout << std::setprecision(3) << "sig2o (var) " << sig2o << " " << sig2o_mean_var << " " << std::endl;
 
    // finalize score_mean
    double sig2o_tilde = sig2o+std::sqrt(sig2o_mean_var);
    score_mean = score_mean / sig2o_tilde;
    // score sample variance
    score_var  = M2a / (p-1) * 1/(sig2o_tilde*sig2o_tilde);
    // score variance of mean (see e.g. Taylor "an Intro to error analysis..." 1982, p. 102)
    score_mean_var = score_var / p;
    if(debug) std::cout << std::setprecision(3) << "score_mean (var) " << score_mean << " " << score_mean_var << " " << std::endl;
    // normal exit
    return true;
}
