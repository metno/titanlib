#include <iostream>
#include <vector>
#include "titanlib.h"
#include <boost/math/distributions/gamma.hpp>

float titanlib::calc_gamma(float shape, float scale) {
    boost::math::gamma_distribution<> dist(shape, scale);
    float value = boost::math::quantile(dist, 0.5);
    return value;
}

titanlib::Dataset::Dataset(fvec ilats, fvec ilons, fvec ielevs, fvec ivalues) {
    lats = ilats;
    lons = ilons;
    elevs = ielevs;
    values = ivalues;
    flags.resize(lats.size());
}

bool titanlib::Dataset::range_check(const fvec min, const fvec max, const ivec indices) {
    if(indices.size() > 0) {
        ivec iflags;
        titanlib::range_check(subset(values, indices), min, max, iflags);
        unsubset(iflags, flags, indices);
    }
    else {
        titanlib::range_check(values, min, max, flags);
    }
}

bool titanlib::Dataset::range_check_climatology(int unixtime, const fvec plus, const fvec minus, const ivec indices) {
    if(indices.size() > 0) {
        ivec iflags;
        titanlib::range_check_climatology(subset(lats, indices), subset(lons, indices), subset(elevs, indices), subset(values, indices), unixtime, plus, minus, iflags);
        unsubset(iflags, flags, indices);
    }
    else {
        titanlib::range_check_climatology(lats, lons, elevs, values, unixtime, plus, minus, flags);
    }
}
bool titanlib::Dataset::sct(int nmin, int nmax, int nminprof, float dzmin, float dhmin, float dz, const fvec t2pos, const fvec t2neg, const fvec eps2, fvec& sct, const ivec indices) {
    if(indices.size() > 0) {
        ivec iflags;
        titanlib::sct(subset(lats, indices), subset(lons, indices), subset(elevs, indices), subset(values, indices), nmin, nmax, nminprof, dzmin, dhmin , dz, subset(t2pos, indices), subset(t2neg, indices), subset(eps2, indices), sct, iflags);
        unsubset(iflags, flags, indices);
    }
    else {
        titanlib::sct(lats, lons, elevs, values, nmin, nmax, nminprof, dzmin, dhmin , dz, t2pos, t2neg, eps2, sct, flags);
    }
}
fvec titanlib::Dataset::subset(const fvec& array, const ivec& indices) {
    fvec new_array = fvec(indices.size());
    for(int i = 0; i < indices.size(); i++) {
        new_array[i] = array[indices[i]];
    }
    return new_array;
}
void titanlib::Dataset::unsubset(const ivec& array, ivec& orig_array, const ivec& indices) {
    assert(array.size() == indices.size());
    for(int i = 0; i < indices.size(); i++) {
        orig_array[indices[i]] = array[i];
    }
}
