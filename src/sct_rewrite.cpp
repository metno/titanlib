#include <vector>
#include <math.h>
#include "titanlib.h"
#include <assert.h>
#include <iostream>
#include <numeric>
#include <gsl/gsl_multimin.h>


// helpers
float compute_quantile(double quantile, const fvec& array);

// forward declarations
fvec compute_vertical_profile(const fvec& lats, const fvec& lons, const fvec& elevs, const fvec& values, double meanT, double gamma, double a, double exact_p10, double exact_p90, int nminprof, double dzmin);
double basic_vertical_profile_optimizer_function(const gsl_vector *v, void *data);
fvec basic_vertical_profile(const int n, const double *elevs, const double t0, const double gamma);
double vertical_profile_optimizer_function(const gsl_vector *v, void *data);
fvec vertical_profile(const int n, const double *elevs, const double t0, const double gamma, const double a, const double h0, const double h1i);


// start SCT //
int spatial_consistency_test(const fvec& lats, const fvec& lons, const fvec& elevs, const fvec& values, ivec& flags, int nminprof, double dzmin)
{
    const int s = values.size();
    // assert that the arrays we expect are of size s
    if( lats.size() != s || lons.size() != s || elevs.size() != s || values.size() != s) { return 1; }

    flags.resize(s, 0);

    /*
    Stuff for VP
    */
    double gamma = -0.0065;
    double a = 5.0;
    double meanT = std::accumulate(values.begin(), values.end(), 0.0) / s;
    double exact_p10 = compute_quantile(0.10, elevs);
    double exact_p90 = compute_quantile(0.90, elevs);

    // calculate background
    fvec vp = compute_vertical_profile(lats, lons, elevs, values, meanT, gamma, a, exact_p10, exact_p90, nminprof, dzmin);

    // now have temperature profile (vp)
    for(int i=0; i<s; i++) {
        assert(vp[i] !=-999);
    }

    // TODO: more SCT stuff to come...


    return 0;
}
// end SCT //


fvec compute_vertical_profile(const fvec& lats, const fvec& lons, const fvec& elevs, const fvec& values, double meanT, double gamma, double a, double exact_p10, double exact_p90, int nminprof, double dzmin) {

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
    std::vector<double> delevs(elevs.begin(), elevs.end());
    std::vector<double> dvalues(values.begin(), values.end());
    double * dpelevs = delevs.data();
    double * dpvalues = dvalues.data();
    double * data[3] = {&nd, dpelevs, dpvalues};

    // Check if terrain is too flat
    double z05 = compute_quantile(0.05, elevs);
    double z95 = compute_quantile(0.95, elevs);

    // should we use the basic or more complicated vertical profile?
    bool use_basic = elevs.size() < nminprof || (z95 - z05) < dzmin;

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
    fvec vp;
    if(use_basic) {
        vp = basic_vertical_profile(nd, dpelevs, gsl_vector_get(s->x, 0), gsl_vector_get(s->x, 1));
    }
    else {
        // then actually calculate the vertical profile using the minima
        vp = vertical_profile(nd, dpelevs,  gsl_vector_get(s->x, 0), gsl_vector_get(s->x, 1),
            gsl_vector_get(s->x, 2), gsl_vector_get(s->x, 3), gsl_vector_get(s->x, 4));
    }

    gsl_vector_free(input);
    gsl_vector_free(ss);
    gsl_multimin_fminimizer_free (s);
    
    return vp;
}

double basic_vertical_profile_optimizer_function(const gsl_vector *v, void *data) 
{
    double **p = (double **)data;
    int n = (int) *p[0]; // is of type double but should be an int
    double *dpelevs = p[1];
    double *dpvalues = p[2];

    // the parameters to mess with
    double meanT = gsl_vector_get(v,0);
    double gamma = gsl_vector_get(v,1);

    // give everything to vp to compute t_out
    fvec t_out = basic_vertical_profile(n, dpelevs, meanT, gamma);
    // RMS
    float total = 0;
    for(int i=0; i<n; i++) {
        total += pow((t_out[i]-dpvalues[i]),2);
    }
    double value = log(pow((total / n),0.5));
    return value;
}

fvec basic_vertical_profile(const int n, const double *elevs, const double t0, const double gamma)
{
    fvec t_out;
    t_out.resize(n, -999);

    for(int i=0; i<n; i++) {
        t_out[i] = t0 + gamma*elevs[i];
    }
    return t_out;
}

double vertical_profile_optimizer_function(const gsl_vector *v, void *data)
{
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
    fvec t_out = vertical_profile(n, dpelevs, meanT, gamma, a, exact_p10, exact_p90);
    // RMS
    double total = 0;
    for(int i=0; i<n; i++) {
        total += pow((t_out[i]-dpvalues[i]),2);
    }
    double value = log(pow((total / n),0.5));
    return value;
}

fvec vertical_profile(const int n, const double *elevs, const double t0, const double gamma, const double a, const double h0, const double h1i)
{
  double h1 = h0 + fabs(h1i); // h1<-h0+abs(h1i)
  // loop over the array of elevations
  fvec t_out;
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
}


//----------------------------------------------------------------------------//
// HELPER FUNCTIONS
float compute_quantile(double quantile, const fvec& array)
{
    int n = array.size();
    fvec array_copy(n);
    // make a copy of the vector
    for(int i = 0; i < n; i++)
        array_copy[i] = array[i];
    float exact_q;
    std::sort(array_copy.begin(), array_copy.end());
    // get the quantile from sorted array
    int lowerIndex = floor(quantile * (n-1));
    int upperIndex = ceil(quantile * (n-1));
    float lowerValue = array_copy[lowerIndex];
    float upperValue = array_copy[upperIndex];
    float lowerQuantile = (float) lowerIndex / (n-1);
    float upperQuantile = (float) upperIndex / (n-1);
    if(lowerIndex == upperIndex) {
        exact_q = lowerValue;
    }
    else {
        assert(upperQuantile > lowerQuantile);
        assert(quantile >= lowerQuantile);
        float f = (quantile - lowerQuantile)/(upperQuantile - lowerQuantile);
        assert(f >= 0);
        assert(f <= 1);
        exact_q = lowerValue + (upperValue - lowerValue) * f;
    }

    return exact_q;
}
