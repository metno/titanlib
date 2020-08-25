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
void titanlib::Dataset::sct( int num_min, 
            int num_max, 
            float inner_radius, 
            float outer_radius,
            int num_iterations,
            int num_min_prof,
            float min_elev_diff,
            float min_horizontal_scale,
            float max_horizontal_scale,
            int kth_closest_obs_horizontal_scale,
            float vertical_scale, 
            const vec& eps2, 
            const vec& tpos_score, 
            const vec& tneg_score, 
            const vec& t_sod,
            const ivec& obs_to_check,
            const vec& background_values,
            std::string background_elab_type, 
            vec& score, 
            vec& rep, 
            vec& sod, 
            vec& num_inner, 
            vec& horizontal_scale, 
            vec& an_inc, 
            vec& an_res, 
            vec& cv_res, 
            vec& innov, 
            vec& idi, 
            vec& idiv, 
            vec& sig2o, 
            const ivec& indices) {
    ivec boxids;
    if(indices.size() > 0) {
        ivec iflags = titanlib::sct(subset(lats, indices), subset(lons, indices), subset(elevs, indices), subset(values, indices), subset(obs_to_check, indices), subset(background_values, indices), background_elab_type, num_min, num_max, inner_radius, outer_radius, num_iterations, num_min_prof, min_elev_diff, min_horizontal_scale , max_horizontal_scale, kth_closest_obs_horizontal_scale, vertical_scale, subset(eps2, indices), subset(tpos_score, indices), subset(tneg_score, indices), subset(t_sod, indices), score, rep, sod, num_inner, horizontal_scale, an_inc, an_res, cv_res, innov, idi, idiv, sig2o);
        unsubset(iflags, flags, indices);
        // DO we have to deal with unsubsetting sct variable?
    } else {
        flags = titanlib::sct( lats, lons, elevs, values, obs_to_check, background_values, background_elab_type, num_min, num_max, inner_radius, outer_radius, num_iterations, num_min_prof, min_elev_diff, min_horizontal_scale , max_horizontal_scale, kth_closest_obs_horizontal_scale, vertical_scale, eps2, tpos_score, tneg_score, t_sod, score, rep, sod, num_inner, horizontal_scale, an_inc, an_res, cv_res, innov, idi, idiv, sig2o);
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
