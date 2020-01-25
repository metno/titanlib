#include <iostream>
#include <vector>
#include "titanlib.h"
#include <boost/math/distributions/gamma.hpp>

float calc_gamma(float shape, float scale) {
    boost::math::gamma_distribution<> dist(shape, scale);
    float value = boost::math::quantile(dist, 0.5);
    return value;
}
bool sct(const vec lats,
        const vec lons,
        const vec elevs,
        const vec values,
        int nmin,
        int nmax,
        int nminprof,
        float dzmin,
        float dhmin,
        float dz,
        const vec t2pos,
        const vec t2neg,
        const vec eps2,
        vec& flags) {
    flags.resize(lats.size());
    for(int i = 0; i < lats.size(); i++) {
        flags[i] = values[i] > 0;
    }
    return true;
}

bool plausibility(const vec lats,
        const vec lons,
        const vec elevs,
        int unixtime,
        float min,
        float max,
        vec& flags) {

}
