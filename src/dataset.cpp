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

//+ SCT - Spatial Consistency Test
void titanlib::Dataset::sct( int num_min_outer, int num_max_outer, float inner_radius, float outer_radius, int num_iterations, int num_min_prof, float min_elev_diff, float min_horizontal_scale, float max_horizontal_scale, int kth_closest_obs_horizontal_scale, float vertical_scale, float value_minp, float value_maxp, const vec& value_mina, const vec& value_maxa, const vec& value_minv, const vec& value_maxv, const vec& eps2, const vec& tpos, const vec& tneg, bool debug, const ivec& obs_to_check, const vec& background_values, std::string background_elab_type, vec& scores, const ivec& indices) {
    ivec boxids;
    if(indices.size() > 0) {
        ivec iflags = titanlib::sct(subset(lats, indices), subset(lons, indices), subset(elevs, indices), subset(values, indices), subset(obs_to_check, indices), subset(background_values, indices), background_elab_type, num_min_outer, num_max_outer, inner_radius, outer_radius, num_iterations, num_min_prof, min_elev_diff, min_horizontal_scale , max_horizontal_scale, kth_closest_obs_horizontal_scale, vertical_scale, value_minp, value_maxp, subset(value_mina, indices), subset(value_maxa, indices), subset(value_minv, indices), subset(value_maxv, indices), subset(eps2, indices), subset(tpos, indices), subset(tneg, indices), debug, scores);
        unsubset(iflags, flags, indices);
        // DO we have to deal with unsubsetting sct variable?
    } else {
        flags = titanlib::sct( lats, lons, elevs, values, obs_to_check, background_values, background_elab_type, num_min_outer, num_max_outer, inner_radius, outer_radius, num_iterations, num_min_prof, min_elev_diff, min_horizontal_scale , max_horizontal_scale, kth_closest_obs_horizontal_scale, vertical_scale, value_minp, value_maxp, value_mina, value_maxa, value_minv, value_maxv, eps2, tpos, tneg, debug, scores);
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
