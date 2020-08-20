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
        // determine if we have too many or too few observations
        // (too many means we can reduce the distance, too few mean isolation problem and cannot flag?)
        int num_min,
        int num_max,
        // first find everything close to the point that we are testing (maxdist)
        float inner_radius,
        float outer_radius,
        int num_iterations,
        int num_min_prof,
        float min_elev_diff,
        float min_horizontal_scale,
        float vertical_scale,
        const vec& pos,
        const vec& neg,
        const vec& eps2,
        vec& prob_gross_error,
        vec& rep) {
/*
 Spatial Consistency Test. Flag observations that are (likely) affected by
  gross measurement errors (or large representativeness errors) based on 
  neighbouring observations.

 Reference: Lussana, C., Uboldi, F. and Salvati, M.R. (2010), A spatial 
             consistency test for surface observations from mesoscale 
             meteorological networks. Q.J.R. Meteorol. Soc., 
             136: 1075-1088. doi:10.1002/qj.622
 (Reference Abbreviation: LUS10)

 -- Notes and Definitions --

  o Loop over s observations (observation index is "curr")
  | + select the subset of s_box closest observations to the curr-th
  |    then perform SCT considering only this subset
  | + observations (s_box-vector): 
  |    (observed_)values = true_values + observation_errors
  |    NOTE: true_values are unknown, useful for theory
  | + background (s_box-vector):
  |    vp = true_value + background_errors
  |    ("vp" stands for vertical profile, because we used sct mostly for temperature)
  | + observation_errors (s_box-vector): 
  |   random variable that follows a multivariate normal (MVN) distribution with:
  |     R = s_box x s_box covariance matrix, R = sig2o * identity_matrix 
  |     mean (s_box-vector) = 0 (for all s elements)
  | + background_errors (s_box-vector):
  |   random variable that follows a MVN distribution with:
  |     S= s_box x s_box covariance matrix, S = exponential 
  |     mean (s_box-vector) = 0
  |
  |
  | A few more definitions:
  | + innovation (s_box-vector)= observations - background = values - vp
  | + analysis residuals (s_box-vector)= observations - analysis 
  | + analysis increments (s_box-vector)= analysis - background
  | + probability of gross error (pog, s_box-vector)
  | + coefficient of representativeness (corep, s_box-vector)
  o

*/
    const int s = values.size();
    if( lats.size() != s || lons.size() != s || elevs.size() != s || values.size() != s || pos.size() != s || neg.size() != s || eps2.size() != s)
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
    if(vertical_scale <= 0)
        throw std::invalid_argument("vertical_scale must be > 0");
    if(inner_radius < 0)
        throw std::invalid_argument("inner_radius must be >= 0");
    if(outer_radius < inner_radius)
        throw std::invalid_argument("outer_radius must be >= inner_radius");

    ivec flags(s, 0);
    prob_gross_error.clear();
    prob_gross_error.resize(s, 0);
    rep.clear();
    rep.resize(s, 0);

    titanlib::KDTree tree(lats, lons);

    for(int iteration = 0; iteration < num_iterations; iteration++) {
        double s_time0 = titanlib::util::clock();

        int thrown_out = 0; // reset this number each loop (this is for breaking if we don't throw anything new out)

        ivec checked(s, 0);  // Keep track of which observations have been checked
        int count_oi = 0;
        for(int curr=0; curr < s; curr++) {
            // break out if station already flagged
            if(flags[curr] != 0) {
                checked[curr] = 1;
                continue;
            }
            if(checked[curr] > 0) {
                continue;
            }
            // get all neighbours that are close enough
            vec distances;
            ivec neighbour_indices = tree.get_neighbours_with_distance(lats[curr], lons[curr], outer_radius, num_max, true, distances);
            neighbour_indices = remove_flagged(neighbour_indices, flags);
            if(neighbour_indices.size() < num_min) {
                checked[curr] = 1;
                // flag as isolated? 
                continue; // go to next station, skip this one
            }

            // call SCT with this box 
            vec lons_box = titanlib::util::subset(lons, neighbour_indices);
            vec elevs_box = titanlib::util::subset(elevs, neighbour_indices);
            vec lats_box = titanlib::util::subset(lats, neighbour_indices);
            vec values_box = titanlib::util::subset(values, neighbour_indices);
            vec eps2_box = titanlib::util::subset(eps2, neighbour_indices);
            int s_box = neighbour_indices.size();
            // the thing to flag is at "curr", ano not included in the box

            // Compute the background
            vec vp;
            if(num_min_prof >= 0) {
                vp = compute_vertical_profile(elevs_box, elevs_box, values_box, num_min_prof, min_elev_diff);
            }
            else {
                double meanT = std::accumulate(values_box.begin(), values_box.end(), 0.0) / values_box.size();
                vp.resize(s, -999);
                for(int l = 0; l < s; l++) {
                    vp[l] = meanT;
                }
            }
           
            /* Compute Dh. 
             The location-dependent horizontal de-correlation lenght scale 
             used for the background error correlation matrix */
            boost::numeric::ublas::matrix<float> disth(s, s);
            boost::numeric::ublas::matrix<float> distz(s, s);
            boost::numeric::ublas::vector<float> Dh(s);

            for(int i=0; i < s_box; i++) {
                vec Dh_vector(s_box);
                for(int j=0; j < s_box; j++) {
                    disth(i, j) = titanlib::util::calc_distance(lats_box[i], lons_box[i], lats_box[j], lons_box[j]);
                    distz(i, j) = fabs(elevs_box[i] - elevs_box[j]);
                    if(i != j) {
                        if(i < j)
                            Dh_vector[j - 1] = disth(i, j);
                        else if(i > j)
                            Dh_vector[j] = disth(i, j);
                    }
                }
                Dh(i) = titanlib::util::compute_quantile(0.10, Dh_vector);
            }

            double Dh_mean = std::accumulate(std::begin(Dh), std::end(Dh), 0.0) / Dh.size();
            if(Dh_mean < min_horizontal_scale) {
                Dh_mean = min_horizontal_scale;
            }
            
            // Compute S + eps2*I and store it in S 
            boost::numeric::ublas::matrix<float> S(s_box,s_box);
            boost::numeric::ublas::matrix<float> Sinv(s_box,s_box);
            for(int i=0; i < s_box; i++) {
                for(int j=0; j < s_box; j++) {
                    double value = std::exp(-.5 * std::pow((disth(i, j) / Dh_mean), 2) - .5 * std::pow((distz(i, j) / vertical_scale), 2));
                    if(i==j) { // weight the diagonal?? (0.5 default)
                        value = value + eps2_box[i];
                    }
                    S(i,j) = value;
                }
            }
            
            // Compute innovations (= observations - background)
            boost::numeric::ublas::vector<float> d(s_box);
            for(int i=0; i < s_box; i++) {
                d(i) = values_box[i] - vp[i]; 
            }

            /* ---------------------------------------------------
               Beginning of real SCT
               ------------------------------------------------------*/
            // Compute ( S + eps2 * I )^(-1) ans store it in Sinv
            bool b = invert_matrix(S, Sinv);
            if(b != true) {
                // TODO: flag differently or give an error???
                continue;
            }

            // Definitions
            boost::numeric::ublas::vector<float> Zinv(s_box), Sinv_d(s_box), ainc(s_box), ares(s_box), cvres(s_box);
            boost::numeric::ublas::vector<float> rep_temp(s_box);
            double sig2o = 0;

            // Matrix multiplications
            for(int i=0; i<s_box; i++) {
                // Recover the actual S from S + eps2*I (unweight the diagonal)
                S(i,i) -= eps2_box[i];
                // Sinv_d = ( S + eps2 * I )^(-1) * innovation
                double acc = 0;
                for(int j=0; j<s_box; j++) {
                    acc += Sinv(i,j)*d(j);
                }
                Sinv_d(i) = acc;
            }
           
            //    analysis increment = analysis - background = S * ( S + eps2 * I )^(-1) * innovation
            //                  Zinv = inverse of the elements on the diagonal of ( S + eps2 * I )^(-1)
            //    analysis residuals = observations - analysis (= innovation - analysis increment)
            // cv-analysis residuals = observations - cv-analysis, LUS10 Eq.(A13)
            // observ error variance = mean( analysis residuals * innovation), LUS10 Eq.(32)
            //                  ppog = quantity proportional to pog
            boost::numeric::ublas::vector<float> ppog(s_box);
            int ccount = 0;
            double M2 = 0;
            double ppog_mean = 0;      // ppog mean
            double ppog_var = 0;       // ppog sample variance
            for(int i=0; i<s_box; i++) {
                double acc = 0;
                for(int j=0; j<s_box; j++) {
                    acc += S(i,j)*Sinv_d(j);
                }
                ainc(i) = acc;
                Zinv(i) = 1/Sinv(i,i); 
                ares(i) = d(i)-ainc(i); 
                cvres(i) = Zinv(i) * Sinv_d(i); 
                rep_temp(i) = d(i)*ares(i);  
                sig2o += rep_temp(i);       
                // Welford's online algorithm for mean and variance 
                float dist = distances[i];
                if(dist <= inner_radius) {
                   ccount++;
                   ppog(i) = cvres(i) * ares(i);
                   float delta = ppog(i) - ppog_mean;
                   ppog_mean += delta / ccount;
                   float delta2 = ppog(i) - ppog_mean;
                   M2 += delta * delta2;
                }
            }
            // finalize variance
            if(ccount >= 2) {
              ppog_var = M2 / (ccount-1);
            } else {
              // one observation in the inner circle, better skip SCT because of possibly large representativeness error?
              continue;
            }
            // ppog variance of mean
            double ppog_mean_var = ppog_var / ccount;

            sig2o = sig2o/s_box;
            if(sig2o < 0.01) { // too small sig2o values are not allowed (shoudl this be a parameter passed to SCT?)
                sig2o = 0.01;
            }

            /* pog = cv-analysis residuals * analysis residuals / sig2o, LUS10 Eq.(22)
               rep = innovation * analysis residuals / sig2o, LUS10 Eq.(32) numerator
               i-th observation is flagged if ( a) i-th pog is larger than a pre-set threshold) AND
               ( b) i-th ppog is an "outlier" compared to the statistics of ppog in the inner circle)
               condition b) helps to avoid flagging representativeness errors as gross errors,
               at least in those regions where many observations are avaialable
               Outlier = values more than 5 standard deviations from the mean (Lanzante. IJC, 1996)
                then, squared-normalized deviations (snd) is:
                 snd = (deviation from the mean)^2 / (spread + uncertainty in the mean) > 25
                 or, snd = (ppog - ppog_mean)**2   / (ppog_var + ppog_mean_var)         > 25
               */
            ccount = 0;
            for(int i = 0; i < s_box; i++) {
                int index = neighbour_indices[i];
                float dist = distances[i];
                if(dist <= inner_radius) {
                    float pog = ppog(i) / sig2o;
                    float ppog_snd = std::pow( ( ppog(i) - ppog_mean), 2) / (ppog_var + ppog_mean_var);
                    float corep = rep_temp(i) / sig2o;
                    prob_gross_error[index] = std::max(pog, prob_gross_error[index]);
                    rep[index] = std::max(corep, rep[index]);
                    if(((cvres(i) > 0 && pog > pos[index]) || (cvres(i) <= 0 && pog > neg[index])) && ppog_snd>25 ) {
                        flags[index] = 1;
                        thrown_out++;
                    }
                    checked[index] = 1;
                    ccount++;
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
        if(flags[indices[i]] == 0 ) {
            removed.push_back(indices[i]);
        }
    }
    return removed;
}
