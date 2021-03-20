#include <vector>
#include <math.h>
#include "titanlib.h"
#include <assert.h>
#include <iostream>
#include <numeric>
#include <exception>
#include <gsl/gsl_multimin.h>
#include <gsl/gsl_linalg.h>

using namespace titanlib;


/*----------------------------------------------------------------------------
 BACKGROUND FUNCTIONS */

namespace {
    vec compute_vertical_profile(const vec& elevs, const vec& oelevs, const vec& values, int num_min_prof, double min_elev_diff, bool debug);

    vec compute_vertical_profile_Theil_Sen(const vec& elevs, const vec& oelevs, const vec& values, int num_min_prof, double min_elev_diff, bool debug);

    double basic_vertical_profile_optimizer_function(const gsl_vector *v, void *data);

    vec basic_vertical_profile(const int n, const double *elevs, const double t0, const double gamma);

    double vertical_profile_optimizer_function(const gsl_vector *v, void *data);

    vec vertical_profile(const int n, const double *elevs, const double t0, const double gamma, const double a, const double h0, const double h1i);
}

/*----------------------------------------------------------------------------
 BACKGROUND FUNCTIONS */

//+ Compute the background for all the observations in the outer circle
vec titanlib::background(const vec& elevs, const vec& values, int num_min_prof, float min_elev_diff, float value_minp, float value_maxp, BackgroundType background_type, const vec& external_background_values, const ivec& indices_global_outer, bool debug) {

    vec bvalues;
    int p = values.size();

    if( background_type == titanlib::VerticalProfile) {
        // vertical profile is independent from the curr-th observation
        bvalues = ::compute_vertical_profile(elevs, elevs, values, num_min_prof, min_elev_diff, debug);
    } else if( background_type == titanlib::VerticalProfileTheilSen) {
        // vertical profile is independent from the curr-th observation
        bvalues = ::compute_vertical_profile_Theil_Sen(elevs, elevs, values, num_min_prof, min_elev_diff, debug);
    } else if( background_type == titanlib::MeanOuterCircle){
        // vertical profile is independent from the curr-th observation
        float mean_val = std::accumulate(values.begin(), values.end(), 0.0) / p;
        bvalues.resize(p, mean_val);
    } else if( background_type == titanlib::MedianOuterCircle){
        // vertical profile is independent from the curr-th observation
        float median_val = titanlib::compute_quantile( 0.5, values);
        bvalues.resize(p, median_val);
    } else if( background_type == titanlib::External){
        bvalues = titanlib::subset(external_background_values, indices_global_outer);
    }
    // set the range 
    for(int i=0; i<p; i++) {
        if(bvalues[i] < value_minp) bvalues[i] = value_minp;
        if(bvalues[i] > value_maxp) bvalues[i] = value_maxp;
    }
    //
    return bvalues;
}
namespace {
vec compute_vertical_profile_Theil_Sen(const vec& elevs, const vec& oelevs, const vec& values, int num_min_prof, double min_elev_diff, bool debug) {
    // Starting value guesses
    double gamma = -0.0065;
    double min_diff = (double) min_elev_diff;
    double meanT = std::accumulate(values.begin(), values.end(), 0.0) / values.size();

    // Special case when all observations have the same elevation
    if ( std::min_element(elevs.begin(),elevs.end()) == std::max_element(elevs.begin(),elevs.end())) {
        vec vp( oelevs.size(), meanT);
        return vp;
    }

    // Check if terrain is too flat
    double z05 = titanlib::compute_quantile(0.05, elevs);
    double z95 = titanlib::compute_quantile(0.95, elevs);

    // should we use the basic or more complicated vertical profile?
    bool use_basic = elevs.size() < num_min_prof || (z95 - z05) < min_diff;

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
        m_median = titanlib::compute_quantile( 0.5, m);
    }
    for(int i=0; i<n; i++) {
        q[i] = values[i] - m_median * elevs[i];
    }
    double q_median = titanlib::compute_quantile( 0.5, q);
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
    double min_diff = (double) min_elev_diff;
    double a = 5.0;

    double meanT = std::accumulate(values.begin(), values.end(), 0.0) / values.size();

    // Special case when all observations have the same elevation
    if ( std::min_element(elevs.begin(),elevs.end()) == std::max_element(elevs.begin(),elevs.end())) {
        vec vp( oelevs.size(), meanT);
        return vp;
    }
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
    bool use_basic = elevs.size() < num_min_prof || (z95 - z05) < min_diff;

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

}
