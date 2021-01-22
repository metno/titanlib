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
    flags.resize(lats.size(), 0);
}

void titanlib::Dataset::range_check(const vec& min, const vec& max, const ivec& indices) {
    ivec new_flags = titanlib::range_check(subset(values, indices), subset(min, indices), subset(max, indices));
    merge(new_flags, subset(indices));
}

void titanlib::Dataset::range_check_climatology(int unixtime, const vec& pos, const vec& neg, const ivec& indices) {
    ivec new_flags = titanlib::range_check_climatology(subset(points, indices), subset(values, indices), unixtime, pos, neg);
    merge(new_flags, subset(indices));
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
    ivec new_flags = titanlib::sct(subset(points, indices), subset(values, indices), num_min, num_max, inner_radius, outer_radius, num_iterations, num_min_prof, min_elev_diff, min_horizontal_scale , vertical_scale, subset(t2pos, indices), subset(t2neg, indices), subset(eps2, indices), prob_gross_error, rep);
    // TODO: Do we have to deal with unsubsetting sct variable?
    merge(new_flags, subset(indices));
}

void titanlib::Dataset::buddy_check(const vec& radius, const ivec& num_min, float threshold, float max_elev_diff, float elev_gradient, float min_std, int num_iterations, const ivec& obs_to_check, const ivec& indices) {
    ivec new_flags = titanlib::buddy_check(subset(points, indices), subset(values, indices), subset(radius, indices), subset(num_min, indices), threshold, max_elev_diff, elev_gradient, min_std, num_iterations, subset(obs_to_check, indices));
    merge(new_flags, subset(indices));
}
void titanlib::Dataset::buddy_event_check(const vec& radius, const ivec& num_min, float event_threshold, float threshold, float max_elev_diff, float elev_gradient, int num_iterations, const ivec& obs_to_check, const ivec& indices) {
    ivec new_flags = titanlib::buddy_event_check(subset(points, indices), subset(values, indices), subset(radius, indices), subset(num_min, indices), event_threshold, threshold, max_elev_diff, elev_gradient, num_iterations, subset(obs_to_check, indices));
    merge(new_flags, subset(indices));
}
void titanlib::Dataset::isolation_check(int num_min, float radius, float vertical_radius, const ivec& indices) {
    ivec new_flags = titanlib::isolation_check(subset(points, indices), num_min, radius, vertical_radius);
    merge(new_flags, subset(indices));
}
void titanlib::Dataset::duplicate_check(float radius, float vertical_range, const ivec& indices) {
    ivec new_flags = titanlib::duplicate_check(subset(points, indices), radius, vertical_range);
    merge(new_flags, subset(indices));
}
void titanlib::Dataset::dem_check(const vec& dem, float max_elev_diff) {
    const vec& elevs = points.get_elevs();
    for(int i = 0; i < points.size(); i++) {
        float elev_diff = fabs(elevs[i] - dem[i]);
        if(elev_diff > max_elev_diff) {
            flags[i] = 1;
        }
    }
}
void titanlib::Dataset::external_check(const ivec& flags) {
    if(flags.size() != this->flags.size())
        throw std::invalid_argument("External flags different size than dataset flags");

    for(int i = 0; i < this->flags.size(); i++) {
        if(flags[i] > 0) {
            this->flags[i] = 1;
        }
    }
}
void titanlib::Dataset::metadata_check(bool check_lat, bool check_lon, bool check_elev, bool check_laf, const ivec& indices) {
    ivec new_flags = titanlib::metadata_check(subset(points, indices), check_lat, check_lon, check_elev, check_laf);
    merge(new_flags, subset(indices));
}
void titanlib::Dataset::merge(const ivec& new_flags, const ivec& indices) {
    if(new_flags.size() != indices.size())
        throw std::invalid_argument("new_flags and indices must be the same length");

    const int s = new_flags.size();
    for(int i = 0; i < s; i++) {
        if(new_flags[i] > 0) {
            int index = indices[i];
            // std::cout << i << " " << index << " " << flags[i] << std::endl;
            if(index >= flags.size() || index < 0)
                throw std::runtime_error("One or more indices are invalid");
            flags[index] = new_flags[i];
        }
    }
}
vec titanlib::Dataset::subset(const vec& array, const ivec& indices, ivec& new_indices) {
    // TODO: Remove this function
    // if(array.size() == 1 && flags.size() != 1)
    //     return array
    vec results;
    results.reserve(array.size());
    new_indices.clear();
    new_indices.reserve(array.size());
    ivec indices0 = indices;
    if(indices0.size() == 0) {
        indices0.resize(array.size());
        for(int i = 0; i < array.size(); i++)
            indices0[i] = i;
    }
    for(int i = 0; i < indices0.size(); i++) {
        int index = indices0[i];
        if(flags[index] == 0) {
            results.push_back(array[index]);
            new_indices.push_back(index);
        }
    }
    return results;
}

ivec titanlib::Dataset::subset(const ivec& indices) {
    ivec indices0 = indices;
    if(indices0.size() == 0) {
        indices0.resize(flags.size());
        for(int i = 0; i < flags.size(); i++)
            indices0[i] = i;
    }
    ivec results;
    results.reserve(indices0.size());
    for(int i = 0; i < indices0.size(); i++) {
        int index = indices0[i];
        if(flags[index] == 0) {
            results.push_back(index);
        }
    }
    return results;
}

vec titanlib::Dataset::subset(const vec& array, const ivec& indices) {
    if(array.size() == 1 && flags.size() != 1)
        return array;
    vec results;
    results.reserve(array.size());
    ivec indices0 = indices;
    if(indices0.size() == 0) {
        indices0.resize(array.size());
        for(int i = 0; i < array.size(); i++)
            indices0[i] = i;
    }
    for(int i = 0; i < indices0.size(); i++) {
        int index = indices0[i];
        if(flags[index] == 0) {
            results.push_back(array[index]);
        }
    }
    return results;
}

Points titanlib::Dataset::subset(const Points& input, const ivec& indices) {
    ivec indices0 = indices;
    if(indices0.size() == 0) {
        indices0.resize(points.size());
        for(int i = 0; i < points.size(); i++)
            indices0[i] = i;
    }

    int size = indices0.size();
    vec ilats = input.get_lats();
    vec ilons = input.get_lons();
    vec ielevs = input.get_elevs();
    vec ilafs = input.get_lafs();
    vec lats;
    vec lons;
    vec elevs;
    vec lafs;
    lats.reserve(size);
    lons.reserve(size);
    elevs.reserve(size);
    lafs.reserve(size);

    for(int i=0; i < size; i++) {
        int index = indices0[i];
        assert(index < flags.size());
        if(flags[index] == 0) {
            lats.push_back(ilats[index]);
            lons.push_back(ilons[index]);
            elevs.push_back(ielevs[index]);
            lafs.push_back(ilafs[index]);
        }
    }
    return Points(lats, lons, elevs, lafs);
}
