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

void titanlib::Dataset::range_check(const vec& min, const vec& max, const ivec& indices) {
    if(indices.size() > 0) {
        ivec iflags = titanlib::range_check(subset(values, indices), subset(min, indices), subset(max, indices));
        unsubset(iflags, flags, indices);
    }
    else {
        flags = titanlib::range_check(values, min, max);
    }
}

void titanlib::Dataset::range_check_climatology(int unixtime, const vec& pos, const vec& neg, const ivec& indices) {
    if(indices.size() > 0) {
        ivec iflags = titanlib::range_check_climatology(subset(lats, indices), subset(lons, indices), subset(elevs, indices), subset(values, indices), unixtime, subset(pos, indices), subset(neg, indices));
        unsubset(iflags, flags, indices);
    }
    else {
        flags = titanlib::range_check_climatology(lats, lons, elevs, values, unixtime, pos, neg);
    }
}

void titanlib::Dataset::buddy_check(const vec& radius, const ivec& num_min, float threshold, float max_elev_diff, float elev_gradient, float min_std, int num_iterations, const ivec& obs_to_check, const ivec& indices) {
    if(indices.size() > 0) {
        ivec iflags = titanlib::buddy_check(lats, lons, elevs, values, subset(radius, indices), subset(num_min, indices), threshold, max_elev_diff, elev_gradient, min_std, num_iterations, obs_to_check);
        unsubset(iflags, flags, indices);
    }
    else {
        flags = titanlib::buddy_check(lats, lons, elevs, values, radius, num_min, threshold, max_elev_diff, elev_gradient, min_std, num_iterations, obs_to_check);
    }
}
void titanlib::Dataset::buddy_event_check(const vec& radius, const ivec& num_min, float event_threshold, float threshold, float max_elev_diff, float elev_gradient, int num_iterations, const ivec& obs_to_check, const ivec& indices) {
    if(indices.size() > 0) {
        ivec iflags = titanlib::buddy_event_check(lats, lons, elevs, values, subset(radius, indices), subset(num_min, indices), event_threshold, threshold, max_elev_diff, elev_gradient, num_iterations, obs_to_check);
        unsubset(iflags, flags, indices);
    }
    else {
        flags = titanlib::buddy_event_check(lats, lons, elevs, values, subset(radius, indices), subset(num_min, indices), event_threshold, threshold, max_elev_diff, elev_gradient, num_iterations, obs_to_check);
    }
}
void titanlib::Dataset::isolation_check(int num_min, float radius, float vertical_radius) {
    flags = titanlib::isolation_check(lats, lons, elevs, num_min, radius, vertical_radius);
}
