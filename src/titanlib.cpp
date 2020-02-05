#include <iostream>
#include <vector>
#include "titanlib.h"
#include <boost/math/distributions/gamma.hpp>

float titanlib::calc_gamma(float shape, float scale) {
    boost::math::gamma_distribution<> dist(shape, scale);
    float value = boost::math::quantile(dist, 0.5);
    return value;
}

bool titanlib::sct(const fvec lats,
        const fvec lons,
        const fvec elevs,
        const fvec values,
        int nmin,
        int nmax,
        int nminprof,
        float dzmin,
        float dhmin,
        float dz,
        const fvec t2pos,
        const fvec t2neg,
        const fvec eps2,
        ivec& flags) {
    flags.resize(lats.size());
    for(int i = 0; i < lats.size(); i++) {
        flags[i] = values[i] > 0;
    }
    return true;
}
