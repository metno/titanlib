#include <iostream>
#include <vector>
#include "titanlib.h"

using namespace titanlib;

titanlib::Dataset::Dataset(Points points, vec ivalues) {
    this->points = points;
    lats = points.get_lats();
    lons = points.get_lons();
    elevs = points.get_elevs();
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
        ivec iflags = titanlib::range_check_climatology(titanlib::subset(points, indices), subset(values, indices), unixtime, subset(pos, indices), subset(neg, indices));
        unsubset(iflags, flags, indices);
    }
    else {
        flags = titanlib::range_check_climatology(points, values, unixtime, pos, neg);
    }
}
void titanlib::Dataset::sct(int num_min, int num_max, float inner_radius, float outer_radius, int num_iterations,
            int num_min_prof,
            float min_elev_diff,
            float min_horizontal_scale,
            float vertical_scale,
            const vec& t2pos,
            const vec& t2neg,
            const vec& eps2,
            vec& prob_gross_error, vec& rep, const ivec& indices) {
    ivec boxids;
    if(indices.size() > 0) {
        ivec iflags = titanlib::sct(titanlib::subset(points, indices), subset(values, indices), num_min, num_max, inner_radius, outer_radius, num_iterations, num_min_prof, min_elev_diff, min_horizontal_scale , vertical_scale, subset(t2pos, indices), subset(t2neg, indices), subset(eps2, indices), prob_gross_error, rep);
        unsubset(iflags, flags, indices);
        // DO we have to deal with unsubsetting sct variable?
    }
    else {
        flags = titanlib::sct(points, values, num_min, num_max, inner_radius, outer_radius, num_iterations, num_min_prof, min_elev_diff, min_horizontal_scale , vertical_scale, t2pos, t2neg, eps2, prob_gross_error, rep);
    }
}

void titanlib::Dataset::buddy_check(const vec& radius, const ivec& num_min, float threshold, float max_elev_diff, float elev_gradient, float min_std, int num_iterations, const ivec& obs_to_check, const ivec& indices) {
    if(indices.size() > 0) {
        ivec iflags = titanlib::buddy_check(points, values, subset(radius, indices), subset(num_min, indices), threshold, max_elev_diff, elev_gradient, min_std, num_iterations, obs_to_check);
        unsubset(iflags, flags, indices);
    }
    else {
        flags = titanlib::buddy_check(points, values, radius, num_min, threshold, max_elev_diff, elev_gradient, min_std, num_iterations, obs_to_check);
    }
}
void titanlib::Dataset::buddy_event_check(const vec& radius, const ivec& num_min, float event_threshold, float threshold, float max_elev_diff, float elev_gradient, int num_iterations, const ivec& obs_to_check, const ivec& indices) {
    if(indices.size() > 0) {
        ivec iflags = titanlib::buddy_event_check(points, values, subset(radius, indices), subset(num_min, indices), event_threshold, threshold, max_elev_diff, elev_gradient, num_iterations, obs_to_check);
        unsubset(iflags, flags, indices);
    }
    else {
        flags = titanlib::buddy_event_check(points, values, subset(radius, indices), subset(num_min, indices), event_threshold, threshold, max_elev_diff, elev_gradient, num_iterations, obs_to_check);
    }
}
void titanlib::Dataset::isolation_check(int num_min, float radius, float vertical_radius) {
    flags = titanlib::isolation_check(points, num_min, radius, vertical_radius);
}
