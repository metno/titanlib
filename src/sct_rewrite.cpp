#include <vector>
#include <math.h>
#include "titanlib.h"
#include <assert.h>
#include <iostream>
#include <boost/math/tools/minima.hpp>


template <class F, class T>
std::pair<T, T> brent_find_minima(F f, T min, T max, int bits);

template <class F, class T>
std::pair<T, T> brent_find_minima(F f, T min, T max, int bits, boost::uintmax_t& max_iter);

// helpers
float compute_quantile(double quantile, const fvec array);

// forward declarations
float basic_vertical_profile_optimizer_function(const fvec elevs, const fvec values, const double meanT);
fvec basic_vertical_profile(const fvec elevs, double t0);
float vertical_profile_optimizer_function(const fvec elevs, const fvec values, const double meanT, const double gamma, const double a, const double exact_p10, const double exact_p90);
fvec vertical_profile(const fvec elevs, const double meanT, const double gamma, const double a, const double h0, const double h1i);

int compute_vertical_profile(const fvec lats, const fvec lons, const fvec elevs,const fvec values, double meanT, double gamma, double a, double exact_p10, double exact_p90, int nminprof, double dzmin, double *vp) {

    // Check if terrain is too flat
    double z05 = compute_quantile(0.05, elevs);
    double z95 = compute_quantile(0.95, elevs);

    // should we use the basic or more complicated vertical profile?
    bool use_basic = elevs.size() < nminprof || (z95 - z05) < dzmin;

    if(use_basic) {
        int bits = std::numeric_limits<double>::digits/2;
        // TODO: unsure how to choose sensible max / min???
        std::pair<float, float> r = brent_find_minima(basic_vertical_profile_optimizer_function(elevs, values, meanT), -4., 4., bits);
    }
    else {

    }

}

float basic_vertical_profile_optimizer_function(const fvec elevs, const fvec values, const double meanT) {

    int n = elevs.size();
    // give everything to vp to compute t_out
    fvec t_out = basic_vertical_profile(elevs, meanT);
    // RMS
    float total = 0;
    for(int i=0; i<n; i++) {
        total += pow((t_out[i]-values[i]),2);
    }
    float value = log(pow((total / n),0.5));
    return value;
}

fvec basic_vertical_profile(const fvec elevs, double t0)
{
    fvec t_out(elevs.size());

    double gamma = -0.0065;
    for(int i=0; i<elevs.size(); i++) {
        t_out[i] = t0 + gamma*elevs[i];
    }
    return t_out;
}

float vertical_profile_optimizer_function(const fvec elevs, const fvec values, const double meanT, const double gamma, const double a, const double exact_p10, const double exact_p90)
{
    int n = elevs.size();
    // give everything to vp to compute t_out
    fvec t_out = vertical_profile(elevs, meanT, gamma, a, exact_p10, exact_p90);
    // RMS
    double total = 0;
    for(int i=0; i<n; i++) {
        total += pow((t_out[i]-values[i]),2);
    }
    double value = log(pow((total / n),0.5));
    return value;
}

fvec vertical_profile(const fvec elevs, const double meanT, const double gamma, const double a, const double h0, const double h1i)
{
  double h1 = h0 + fabs(h1i); // h1<-h0+abs(h1i)
  // loop over the array of elevations (z)
  int n = elevs.size();
  fvec t_out(elevs.size());
  for(int i=0; i<n; i++) {
    // define some bools
    bool z_le_h0 = elevs[i] <= h0; // z.le.h0<-which(z<=h0)
    bool z_ge_h1 = elevs[i] >= h1; // z.ge.h1<-which(z>=h1)
    bool z_in = (elevs[i]>h0 && elevs[i]<h1); // z.in<-which(z>h0 & z<h1)
    if(z_le_h0) {
      t_out[i] = meanT-gamma*elevs[i]-a;
    }
    if(z_ge_h1) {
      t_out[i] = meanT-gamma*elevs[i];
    }
    if(z_in) {
      t_out[i] = meanT-gamma*elevs[i]-a/2*(1+cos(M_PI*(elevs[i]-h0)/(h1-h0)));
    }
  }
}


//----------------------------------------------------------------------------//
// HELPER FUNCTIONS
float compute_quantile(double quantile, const fvec array)
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
