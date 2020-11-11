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

using namespace titanlib;

namespace {
    template<class T1, class T2> struct sort_pair_first {
        bool operator()(const std::pair<T1,T2>&left, const std::pair<T1,T2>&right) {
            return left.first < right.first;
        };
    };
}
// helpers
bool invert_matrix (const boost::numeric::ublas::matrix<float>& input, boost::numeric::ublas::matrix<float>& inverse);
void remove_flagged(ivec indices, vec distances, ivec flags);
vec compute_vertical_profile(const vec& elevs, const vec& oelevs, const vec& values, int num_min_prof, double min_elev_diff);
double basic_vertical_profile_optimizer_function(const gsl_vector *v, void *data);
vec basic_vertical_profile(const int n, const double *elevs, const double t0, const double gamma);
double vertical_profile_optimizer_function(const gsl_vector *v, void *data);
vec vertical_profile(const int n, const double *elevs, const double t0, const double gamma, const double a, const double h0, const double h1i);

// start SCT //
ivec titanlib::sct(const Points& points,
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

    const vec lats = points.get_lats();
    const vec lons = points.get_lons();
    const vec elevs = points.get_elevs();

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
        double s_time0 = titanlib::clock();

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
            ivec neighbour_indices = tree.get_neighbours_with_distance(lats[curr], lons[curr], outer_radius, distances, false);
            remove_flagged(neighbour_indices, distances, flags);

            if(neighbour_indices.size() > num_max) {

                int N = neighbour_indices.size();
                std::vector<std::pair<float,int> > pairs(N);
                for(int i = 0; i < neighbour_indices.size(); i++) {
                    pairs[i] = std::pair<float, int>(distances[i], neighbour_indices[i]);
                }
                std::sort(pairs.begin(), pairs.end(), ::sort_pair_first<float,int>());
                distances.clear();
                neighbour_indices.clear();
                distances.resize(num_max);
                neighbour_indices.resize(num_max);
                for(int i = 0; i < num_max; i++) {
                    distances[i] = pairs[i].first;
                    neighbour_indices[i] = pairs[i].second;
                }
            }
            if(neighbour_indices.size() < num_min) {
                checked[curr] = 1;
                // flag as isolated? 
                continue; // go to next station, skip this one
            }

            // call SCT with this box 
            vec lons_box = titanlib::subset(lons, neighbour_indices);
            vec elevs_box = titanlib::subset(elevs, neighbour_indices);
            vec lats_box = titanlib::subset(lats, neighbour_indices);
            vec values_box = titanlib::subset(values, neighbour_indices);
            vec eps2_box = titanlib::subset(eps2, neighbour_indices);
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

            boost::numeric::ublas::matrix<float> disth(s, s);
            boost::numeric::ublas::matrix<float> distz(s, s);
            boost::numeric::ublas::vector<float> Dh(s);

            for(int i=0; i < s_box; i++) {
                vec Dh_vector(s_box);
                for(int j=0; j < s_box; j++) {
                    disth(i, j) = titanlib::calc_distance(lats_box[i], lons_box[i], lats_box[j], lons_box[j]);
                    distz(i, j) = fabs(elevs_box[i] - elevs_box[j]);
                    if(i != j) {
                        if(i < j)
                            Dh_vector[j - 1] = disth(i, j);
                        else if(i > j)
                            Dh_vector[j] = disth(i, j);
                    }
                }
                Dh(i) = titanlib::compute_quantile(0.10, Dh_vector);
            }

            double Dh_mean = std::accumulate(std::begin(Dh), std::end(Dh), 0.0) / Dh.size();
            if(Dh_mean < min_horizontal_scale) {
                Dh_mean = min_horizontal_scale;
            }

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

            boost::numeric::ublas::vector<float> d(s_box);
            for(int i=0; i < s_box; i++) {
                d(i) = values_box[i] - vp[i]; // difference between actual temp and temperature from vertical profile
            }

            /* ---------------------------------------------------
               Beginning of real SCT
               ------------------------------------------------------*/
            bool b = invert_matrix(S, Sinv);
            if(b != true) {
                // TODO: flag differently or give an error???
                continue;
            }

            // Unweight the diagonal
            for(int i=0; i < s_box; i++) {
                S(i,i) -= eps2_box[i];
            }

            boost::numeric::ublas::vector<float> Zinv(s_box), Sinv_d(s_box), ares_temp(s_box), ares(s_box);
            for(int i=0; i<s_box; i++) {
                double acc = 0;
                for(int j=0; j<s_box; j++) {
                    acc += Sinv(i,j)*d(j);
                }
                Sinv_d(i) = acc;
            }
            for(int i=0; i<s_box; i++) {
                double acc = 0;
                for(int j=0; j<s_box; j++) {
                    acc += S(i,j)*Sinv_d(j);
                }
                ares_temp(i) = acc;
            }
            for(int i=0; i<s_box; i++) {
                Zinv(i) = (1/Sinv(i,i)); //Zinv<-1/diag(SRinv)
                ares(i) = ares_temp(i)-d(i); // ares<-crossprod(S,SRinv.d)-d[sel]
            }

            boost::numeric::ublas::vector<float> cvres(s_box);
            for(int i=0; i<s_box; i++) {
                cvres(i) = -1*Zinv(i) * Sinv_d(i);
            }

            double sig2o = 0;
            boost::numeric::ublas::vector<float> sig2o_temp(s_box), negAres_temp(s_box);
            for(int i=0; i<s_box; i++) {
                negAres_temp(i)=-1*ares(i);
                sig2o_temp(i) = d(i)*negAres_temp(i);
                sig2o += sig2o_temp(i);
            }

            sig2o = sig2o/s_box;
            if(sig2o < 0.01) {
                sig2o = 0.01;
            }

            // boost::numeric::ublas::vector<float> pog(s_box);
            // for(int i=0; i<s_box; i++) {
            //     pog(i) = cvres(i)*ares(i) / sig2o;
            // }
            int ccount = 0;
            for(int i = 0; i < s_box; i++) {
                int index = neighbour_indices[i];
                float dist = distances[i];
                if(dist <= inner_radius) {
                    float pog = cvres(i) * ares(i) / sig2o;
                    prob_gross_error[index] = std::max(pog, prob_gross_error[index]);
                    if((cvres(i) < 0 && pog > pos[index]) || (cvres(i) >= 0 && pog > neg[index])) {
                        flags[index] = 1;
                        thrown_out++;
                    }
                    checked[index] = 1;
                    ccount++;
                }
            }
            count_oi++;
        }
        if(thrown_out == 0) {
            if(iteration + 1 < num_iterations)
                std::cout << "Stopping early after " << iteration + 1<< " iterations" << std::endl;
            break;
        }
        std::cout << "Removing " << thrown_out << " Number of OI " << count_oi << std::endl;
        double e_time0 = titanlib::clock();
        std::cout << e_time0 - s_time0 << std::endl;
    }

    return flags;
}
// end SCT //


vec compute_vertical_profile(const vec& elevs, const vec& oelevs, const vec& values, int num_min_prof, double min_elev_diff) {
    // Starting value guesses
    double gamma = -0.0065;
    double a = 5.0;

    double meanT = std::accumulate(values.begin(), values.end(), 0.0) / values.size();
    double exact_p10 = titanlib::compute_quantile(0.10, elevs);
    double exact_p90 = titanlib::compute_quantile(0.90, elevs);

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
    double z05 = titanlib::compute_quantile(0.05, elevs);
    double z95 = titanlib::compute_quantile(0.95, elevs);

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


void remove_flagged(ivec indices, vec distances, ivec flags) {
    ivec indices_new;
    vec distances_new;
    indices_new.reserve(indices.size());
    distances_new.reserve(indices.size());
    for(int i=0; i<indices.size(); i++) {
        if(flags[indices[i]] == 0 ) {
            indices_new.push_back(indices[i]);
            distances_new.push_back(distances[i]);
        }
    }
}
