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
        const vec& eps2,
        const vec& tpos_score,
        const vec& tneg_score,
        const vec& t_sod,
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
 + Observations (s_box-vector): 
    (observed_)values = true_values + observation_errors
    NOTE: true_values are unknown, useful for theory
 + Background (s_box-vector):
    background_value = true_value + background_errors
 + Observation_errors (s_box-vector): 
   random variable that follows a multivariate normal (MVN) distribution with:
     R = s_box x s_box covariance matrix, R = sig2o * identity_matrix 
     mean (s_box-vector) = 0 (for all s elements)
 + Background_errors (s_box-vector):
   random variable that follows a MVN distribution with:
     S= s_box x s_box covariance matrix, S = exponential 
     mean (s_box-vector) = 0
 + Innovation (s_box-vector)= observations - background 
 + Analysis residuals (s_box-vector)= observations - analysis 
 + Analysis increments (s_box-vector)= analysis - background
 + Cross-validation indicates leave-one-out cross-validation
 + (SCT) Score (score, s_box-vector)
   score = cv-analysis residuals * analysis residuals / sig2o, LUS10 Eq.(22)
 + Coefficient of representativeness (corep, s_box-vector)
   rep = innovation * analysis residuals / sig2o, LUS10 Eq.(32) numerator
 + Spatial Outlier Detection (SOD) score (sod, s_box-vector)
    used to distinguish between acceptable REs and GEs
     sod = (deviation from SCT-score areal average)^2 / spread 

 Algorithm: 
  Loop over the s observations (observation index is "curr"). Select the subset of s_box closest observations to the curr-th observation, then perform SCT considering only this subset.
  curr-th observation is flagged if curr-th score is a) larger than a pre-set threshold AND b) an outlier compared to the statistics of score inside the inner circle
  condition b) outlier if sod>threshold. Helps to avoid flagging REs as GEs, at least in those regions where enough observations are avaialable

 Returned values:
  flags. -9999 = not checked; 0 = passed (good); 1 = failed (bad); 11 = isolated (<2 inside inner); 12 = isolated (<num_min inside outer)

*/
    bool debug = true;

    const int s = values.size();
    if( lats.size() != s || lons.size() != s || elevs.size() != s || values.size() != s || tpos_score.size() != s || tneg_score.size() != s || eps2.size() != s || background_values.size() != s)
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
    if(background_elab_type != "vertical_profile" && background_elab_type != "mean_outer_circle" && background_elab_type != "external")
        throw std::invalid_argument("background_elab_type must be one of vertical_profile, mean_outer_circle or external");
    if(background_elab_type == "vertical_profile" && num_min_prof<0)
        throw std::invalid_argument("num_min_prof must be >=0");
    if( (obs_to_check.size() != s && obs_to_check.size() != 1 && obs_to_check.size() !=0) ) {
        throw std::invalid_argument("'obs_to_check' has an invalid length");
    }
    
    // initializations
    ivec flags(s, -9999.);
    score.clear();
    score.resize(s, -9999.);
    sod.clear();
    sod.resize(s, -9999.);
    num_inner.clear();
    num_inner.resize(s, -9999.);
    horizontal_scale.clear();
    horizontal_scale.resize(s, -9999.);
    an_inc.clear();
    an_inc.resize(s, -9999.);
    an_res.clear();
    an_res.resize(s, -9999.);
    cv_res.clear();
    cv_res.resize(s, -9999.);
    innov.clear();
    innov.resize(s, -9999.);
    idi.clear();
    idi.resize(s, -9999.);
    idiv.clear();
    idiv.resize(s, -9999.);
    sig2o.clear();
    sig2o.resize(s, -9999.);
    rep.clear();
    rep.resize(s, -9999.);

    // if obs_to_check is empty then check all
    bool check_all = obs_to_check.size() != s;

    // KDtree has to do with fast computation of distances
    titanlib::KDTree tree(lats, lons);

    // SCT iterations
    for(int iteration = 0; iteration < num_iterations; iteration++) {
        double s_time0 = titanlib::util::clock();

        int thrown_out = 0; // reset this number each loop (this is for breaking if we don't throw anything new out)

        ivec checked(s, 0);  // Keep track of which observations have been checked
        int count_oi = 0;

        // loop over observations
        for(int curr=0; curr < s; curr++) {
            if(debug) {
              std::cout << "===> curr " << curr << "===============" << std::endl;
            }
            if(!check_all && obs_to_check[curr] != 1) {
                if(debug) {
                  std::cout << "..not to check " << curr << std::endl;
                }
                continue;
            }
            // break out if station already flagged
            if(flags[curr] == 1) {
                checked[curr] = 1;
                if(debug) {
                  std::cout << "..checked1 " << curr << std::endl;
                }
                continue;
            }
            if(checked[curr] > 0) {
                if(debug) {
                    std::cout << "..checked2 " << curr << std::endl;
                }
                continue;
            }
            // no neighbours within the inner circle = distinction between RE and GE not reliable
            vec distances;
            ivec neighbour_indices = tree.get_neighbours_with_distance(lats[curr], lons[curr], inner_radius, 2, true, distances);
            neighbour_indices = remove_flagged(neighbour_indices, flags);
            if(neighbour_indices.size() < 2) {
                flags[curr] = 11;
                checked[curr] = 1;
                if(debug) {
                    std::cout << "@@isolated (inner) " << curr << std::endl;
                }
                continue; // go to next station, skip this one
            }
            // get all neighbours that are close enough
            neighbour_indices = tree.get_neighbours_with_distance(lats[curr], lons[curr], outer_radius, num_max, true, distances);
            neighbour_indices = remove_flagged(neighbour_indices, flags);
            if(neighbour_indices.size() < num_min) {
                flags[curr] = 12;
                checked[curr] = 1;
                if(debug) {
                    std::cout << "@@isolated (outer) " << curr << std::endl;
                }
                continue; // go to next station, skip this one
            }

            // call SCT with this box 
            vec lons_box = titanlib::util::subset(lons, neighbour_indices);
            vec elevs_box = titanlib::util::subset(elevs, neighbour_indices);
            vec lats_box = titanlib::util::subset(lats, neighbour_indices);
            vec values_box = titanlib::util::subset(values, neighbour_indices);
            vec eps2_box = titanlib::util::subset(eps2, neighbour_indices);
            int s_box = neighbour_indices.size();
            if(debug) {
                std::cout << "s_box " << s_box << std::endl;
            }
            // the thing to flag is at "curr", ano not included in the box

            // Compute the background
            vec bvalues_box;
            if( background_elab_type == "vertical_profile") {
                bvalues_box = compute_vertical_profile(elevs_box, elevs_box, values_box, num_min_prof, min_elev_diff);
            } else if( background_elab_type == "mean_outer_circle"){
                double mean_val = std::accumulate(values_box.begin(), values_box.end(), 0.0) / values_box.size();
                bvalues_box.resize(s_box, -9999);
                for(int i = 0; i < s_box; i++) {
                    bvalues_box[i] = mean_val;
                }
            } else if( background_elab_type == "external"){
                bvalues_box = titanlib::util::subset(background_values, neighbour_indices);
            } 
           
            /* Compute Dh. 
             The location-dependent horizontal de-correlation lenght scale 
             used for the background error correlation matrix */
            boost::numeric::ublas::matrix<float> disth(s_box, s_box);
            boost::numeric::ublas::matrix<float> distz(s_box, s_box);
            boost::numeric::ublas::vector<float> Dhbox(s_box);

            for(int i=0; i < s_box; i++) {
                vec Dhbox_vector(s_box);
                for(int j=0; j < s_box; j++) {
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
                Dhbox_mean = min_horizontal_scale;
            }
            if(Dhbox_mean > max_horizontal_scale) {
                Dhbox_mean = max_horizontal_scale;
            }
            if(debug) {
                std::cout << "Dhbox_mean " << Dhbox_mean << std::endl;
            }
            
            // Compute S + eps2*I and store it in S 
            boost::numeric::ublas::matrix<float> S(s_box,s_box);
            boost::numeric::ublas::matrix<float> Sinv(s_box,s_box);
            for(int i=0; i < s_box; i++) {
                for(int j=0; j < s_box; j++) {
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
                if(debug) {
                    std::cout << "ooo Problem in matrix inversion ooo " << curr << std::endl;
                }
                continue;
            }

            // Definitions
            boost::numeric::ublas::vector<float> Zinv(s_box), Sinv_d(s_box), Sinv_1(s_box), ainc(s_box), ares(s_box), cvres(s_box), d(s_box), rep_temp(s_box), idibox(s_box), idivbox(s_box);
            double sig2obox = 0;

            // Matrix multiplications
            // Compute innovations (= observations - background)
            for(int i=0; i<s_box; i++) {
                d(i) = values_box[i] - bvalues_box[i]; 
                // Recover the actual S from S + eps2*I (unweight the diagonal)
                S(i,i) -= eps2_box[i];
                // Sinv_d = ( S + eps2 * I )^(-1) * innovation
                double acc = 0;
                double acc1 = 0;
                for(int j=0; j<s_box; j++) {
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
            boost::numeric::ublas::vector<float> scorebox_numerator(s_box);
            int ccount = 0;
            double M2 = 0;
            double scorebox_mean = 0;      // scorebox mean
            double scorebox_var = 0;       // scorebox sample variance
            for(int i=0; i<s_box; i++) {
                double acc = 0;
                double acc1 = 0;
                for(int j=0; j<s_box; j++) {
                    acc += S(i,j)*Sinv_d(j);
                    acc1 += S(i,j)*Sinv_1(j);
                }
                idibox(i) = acc1;
                ainc(i) = acc;
                Zinv(i) = 1 / Sinv(i,i); 
                ares(i) = d(i) - ainc(i); 
                cvres(i) = Zinv(i) * Sinv_d(i); 
                idivbox(i) = 1 - Zinv(i) * Sinv_1(i); 
                rep_temp(i) = d(i) * ares(i);  
                sig2obox += rep_temp(i);       
                // Welford's online algorithm for mean and variance 
                float dist = distances[i];
                if(dist <= inner_radius) {
                   ccount++;
                   scorebox_numerator(i) = cvres(i) * ares(i);
                   float delta = scorebox_numerator(i) - scorebox_mean;
                   scorebox_mean += delta / ccount;
                   float delta2 = scorebox_numerator(i) - scorebox_mean;
                   M2 += delta * delta2;
                }
            }

            sig2obox = fabs(sig2obox)/s_box;
            if(sig2obox < 0.01) { // negative and too small sig2obox values are not allowed 
                sig2obox = 0.01;
            }
            if(debug) {
                std::cout << "sig2obox " << sig2obox << std::endl;
            }
 
            // finalize scorebox_var and scorebox_mean (and its variance)
            if(ccount >= 2) {
              scorebox_mean = scorebox_mean / sig2obox;
              scorebox_var  = M2 / (ccount-1) * 1/(sig2obox*sig2obox);
            } else {
              // one observation in the inner circle, better skip SCT because of possibly large RE
              // note: the algorithm will never enter this branch, because of check which sets flags=11
              if(debug) {
                  std::cout << "xxx only one obs within inner circle xxx " << std::endl;
              }
              continue;
            }
            // scorebox variance of mean (see e.g. Taylor "an Intro to error analysis..." 1982, p. 102)
            double scorebox_mean_var = scorebox_var / ccount;

            /* score = cv-analysis residuals * analysis residuals / sig2o, LUS10 Eq.(22)
               rep = innovation * analysis residuals / sig2o, LUS10 Eq.(32) numerator
               i-th observation is flagged if i-th score is a) larger than a pre-set threshold AND b) an outlier compared to the statistics of score inside the inner circle
               condition b) helps to avoid flagging REs as GEs, at least in those regions where enough observations are avaialable
               sod (spatial outlier detection score) = (deviation from the mean)^2 / (spread + uncertainty in the mean)
               outlier if sod>threshold 
               */
            for(int i = 0; i < s_box; i++) {
                int index = neighbour_indices[i];
                float dist = distances[i];
                if(dist <= inner_radius) {
                    if(debug) {
                        std::cout << "index " << index << std::endl;
                    }
                    float scorebox = scorebox_numerator(i) / sig2obox;
                    float sodbox = std::pow( ( scorebox_numerator(i) - scorebox_mean), 2) / (scorebox_var + scorebox_mean_var);
                    float corep = rep_temp(i) / sig2obox;
                    // condition defining bad observations 
                    if( ((cvres(i) > 0 && scorebox > tpos_score[index]) || (cvres(i) <= 0 && scorebox > tneg_score[index])) && sodbox > t_sod[index] ) {
                        score[index] = scorebox;
                        rep[index] = corep;
                        sod[index] = sodbox;
                        num_inner[index] = ccount;
                        horizontal_scale[index] = Dhbox_mean;
                        an_inc[index] = ainc(i);
                        an_res[index] = ares(i);
                        cv_res[index] = cvres(i);
                        innov[index] = d(i);
                        idi[index] = idibox(i);
                        idiv[index] = idivbox(i);
                        sig2o[index] = sig2obox;
                        flags[index] = 1;
                        if(debug) {
                            std::cout << std::setprecision(3) << "~~~ index flags innov ares ainc cvres scorebox sod rep: " << index << " " << flags[index] << " " << innov[index] << " " << an_res[index] << " " << an_inc[index] << " " << cv_res[index] << " " << score[index] << " " << sod[index] << " " << rep[index] << " " << std::endl;
                        }
                        thrown_out++;
                    } else {
                        if(scorebox>score[index]) {
                            score[index] = scorebox;
                            rep[index] = corep;
                            sod[index] = sodbox;
                            num_inner[index] = ccount;
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
                                std::cout << std::setprecision(3) << "--- index flags innov ares ainc cvres scorebox sod rep: " << index << " " << flags[index] << " " << innov[index] << " " << an_res[index] << " " << an_inc[index] << " " << cv_res[index] << " " << score[index] << " " << sod[index] << " " << rep[index] << " " << std::endl;
                            }
                        }
                    }
                    checked[index] = 1;
                }
            }
            count_oi++;
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
