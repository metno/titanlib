#include <iostream>
#include <vector>
#include "titanlib.h"

titanlib::Dataset::Dataset(vec ilats, vec ilons, vec ielevs, vec ivalues) {
    lats = ilats;
    lons = ilons;
    elevs = ielevs;
    values = ivalues;
    flags.resize(lats.size());
}

bool titanlib::Dataset::range_check(const vec& min, const vec& max, const ivec indices) {
    int status;
    if(indices.size() > 0) {
        ivec iflags = titanlib::range_check(subset(values, indices), subset(min, indices), subset(max, indices));
        unsubset(iflags, flags, indices);
    }
    else {
        flags = titanlib::range_check(values, min, max);
    }
    if(status == 0) {
        return true;
    }
    else return false;
}

bool titanlib::Dataset::range_check_climatology(int unixtime, const vec& pos, const vec& neg, const ivec indices) {
    int status;
    if(indices.size() > 0) {
        ivec iflags = titanlib::range_check_climatology(subset(lats, indices), subset(lons, indices), subset(elevs, indices), subset(values, indices), unixtime, subset(pos, indices), subset(neg, indices));
        unsubset(iflags, flags, indices);
    }
    else {
        flags = titanlib::range_check_climatology(lats, lons, elevs, values, unixtime, pos, neg);
    }
    if(status == 0) {
        return true;
    }
    else return false;
}
bool titanlib::Dataset::sct(int num_min, int num_max, double inner_radius, double outer_radius, int num_iterations,
            int num_min_prof,
            double dzmin,
            double dhmin,
            float dz,
            const vec& t2pos,
            const vec& t2neg,
            const vec& eps2,
            vec& sct, vec& rep, const ivec indices) {
    int status;
    ivec boxids;
    if(indices.size() > 0) {
        ivec iflags = titanlib::sct(subset(lats, indices), subset(lons, indices), subset(elevs, indices), subset(values, indices), num_min, num_max, inner_radius, outer_radius, num_iterations, num_min_prof, dzmin, dhmin , dz, subset(t2pos, indices), subset(t2neg, indices), subset(eps2, indices), sct, rep);
        unsubset(iflags, flags, indices);
        // DO we have to deal with unsubsetting sct variable?
    }
    else {
        flags = titanlib::sct(lats, lons, elevs, values, num_min, num_max, inner_radius, outer_radius, num_iterations, num_min_prof, dzmin, dhmin , dz, t2pos, t2neg, eps2, sct, rep);
    }
    if(status == 0) {
        return true;
    }
    else return false;
}

bool titanlib::Dataset::buddy_check(const vec& radius, const ivec& buddies_min, const vec& thresholds, float diff_elev_max, float elev_gradient, float min_std, int num_iterations, const ivec& obs_to_check, const ivec indices) {
    int status;
    if(indices.size() > 0) {
        ivec iflags = titanlib::buddy_check(lats, lons, elevs, values, subset(radius, indices), subset(buddies_min, indices), subset(thresholds, indices), diff_elev_max, elev_gradient, min_std, num_iterations, obs_to_check);
        unsubset(iflags, flags, indices);
    }
    else {
        flags = titanlib::buddy_check(lats, lons, elevs, values, radius, buddies_min, thresholds, diff_elev_max, elev_gradient, min_std, num_iterations, obs_to_check);
    }
    if(status == 0) {
        return true;
    }
    else return false;
}
bool titanlib::Dataset::isolation_check(int num_min, float radius, float dz) {
    flags = titanlib::isolation_check(lats, lons, elevs, num_min, radius, dz);
}
