#include <iostream>
#include <vector>
#include "titanlib.h"

titanlib::Dataset::Dataset(fvec ilats, fvec ilons, fvec ielevs, fvec ivalues) {
    lats = ilats;
    lons = ilons;
    elevs = ielevs;
    values = ivalues;
    flags.resize(lats.size());
}

bool titanlib::Dataset::range_check(const fvec& min, const fvec& max, const ivec indices) {
    int status;
    if(indices.size() > 0) {
        ivec iflags;
        status = titanlib::range_check(subset(values, indices), subset(min, indices), subset(max, indices), iflags);
        unsubset(iflags, flags, indices);
    }
    else {
        status = titanlib::range_check(values, min, max, flags);
    }
    if(status == 0) {
        return true;
    }
    else return false;
}

bool titanlib::Dataset::range_check_climatology(int unixtime, const fvec& plus, const fvec& minus, const ivec indices) {
    int status;
    if(indices.size() > 0) {
        ivec iflags;
        status = titanlib::range_check_climatology(subset(lats, indices), subset(lons, indices), subset(elevs, indices), subset(values, indices), unixtime, subset(plus, indices), subset(minus, indices), iflags);
        unsubset(iflags, flags, indices);
    }
    else {
        status = titanlib::range_check_climatology(lats, lons, elevs, values, unixtime, plus, minus, flags);
    }
    if(status == 0) {
        return true;
    }
    else return false;
}
bool titanlib::Dataset::sct(int nmin, int nmax, int nminprof, float dzmin, float dhmin, float dz, const fvec& t2pos, const fvec& t2neg, const fvec& eps2, fvec& sct, fvec& rep, const ivec indices) {
    int status;
    if(indices.size() > 0) {
        ivec iflags;
        status = titanlib::sct(subset(lats, indices), subset(lons, indices), subset(elevs, indices), subset(values, indices), nmin, nmax, nminprof, dzmin, dhmin , dz, subset(t2pos, indices), subset(t2neg, indices), subset(eps2, indices), sct, rep, iflags);
        unsubset(iflags, flags, indices);
        // DO we have to deal with unsubsetting sct variable?
    }
    else {
        status = titanlib::sct(lats, lons, elevs, values, nmin, nmax, nminprof, dzmin, dhmin , dz, t2pos, t2neg, eps2, sct, rep, flags);
    }
    if(status == 0) {
        return true;
    }
    else return false;
}

bool titanlib::Dataset::buddy_check(const fvec& radius, const ivec& buddies_min, const fvec& thresholds, float diff_elev_max, float elev_gradient, float min_std, int num_iterations, const ivec& obs_to_check, const ivec indices) {
    int status;
    if(indices.size() > 0) {
        ivec iflags;
        status = titanlib::buddy_check(lats, lons, elevs, values, subset(radius, indices), subset(buddies_min, indices), subset(thresholds, indices), diff_elev_max, elev_gradient, min_std, num_iterations, flags, obs_to_check);
        unsubset(iflags, flags, indices);
    }
    else {
        status = titanlib::buddy_check(lats, lons, elevs, values, radius, buddies_min, thresholds, diff_elev_max, elev_gradient, min_std, num_iterations, flags, obs_to_check);
    }
    if(status == 0) {
        return true;
    }
    else return false;
}
