#include <iostream>
#include <vector>
#include "titanlib.h"

using namespace titanlib;

titanlib::Dataset::Dataset(Points points, vec ivalues) {
    if(points.size() != ivalues.size()) {
        std::stringstream ss;
        ss << "Points (" << points.size() << ") have a different size than values (" << ivalues.size() <<")";
        throw std::invalid_argument(ss.str());
    }
    this->points = points;
    lats = points.get_lats();
    lons = points.get_lons();
    elevs = points.get_elevs();
    values = ivalues;
    flags.resize(lats.size(), 0);
}

// NOTE: There are two ways to deal with 'indicies', which depends on the kind of quality control
// function:
// Scenario 1: The test for a specific station does not depend on other stations. In this case, we
//             only want to run the test on stations in indices and that are not already flagged.
//             Run subset_valid on all arrays and then use merge_simple to merge flags
// Scenario 2: The test for a specific station does depend on other stations. In this case, we
//             want to run the test on all non-flagged stations (even those not in indices) and then
//             only merge in those specified by indices. Run get_flagged* functions on all arrays, 
//             and use the merge function to merge flags.

void titanlib::Dataset::range_check(const vec& min, const vec& max, const ivec& indices) {
    vec v1 = subset_valid(values, indices);
    vec v2 = subset_valid(min, indices);
    vec v3 = subset_valid(max, indices);
    ivec new_flags = titanlib::range_check(v1, v2, v3);
    // ivec new_flags = titanlib::range_check(subset_valid(values, indices), subset_valid(min, indices), subset_valid(max, indices));
    merge_simple(new_flags, subset_valid(indices));
}

void titanlib::Dataset::range_check_climatology(int unixtime, const vec& pos, const vec& neg, const ivec& indices) {
    ivec new_flags = titanlib::range_check_climatology(subset_valid(points, indices), subset_valid(values, indices), unixtime, pos, neg);
    merge_simple(new_flags, subset_valid(indices));
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
    // Only keep flags that are not flagged
    ivec new_flags = titanlib::sct(get_unflagged_points(), get_unflagged(values), num_min, num_max, inner_radius, outer_radius, num_iterations, num_min_prof, min_elev_diff, min_horizontal_scale , vertical_scale, get_unflagged(t2pos), get_unflagged(t2neg), get_unflagged(eps2), prob_gross_error, rep);
    // ivec new_indices = get_unflagged_indices();
    // ivec new_flags2 = titanlib::subset(new_flags, new_indices);
    // merge_simple(new_flags2, new_indices);
    merge(new_flags, indices);
}

void titanlib::Dataset::buddy_check(const vec& radius, const ivec& num_min, float threshold, float max_elev_diff, float elev_gradient, float min_std, int num_iterations, const ivec& obs_to_check, const ivec& indices) {
    ivec new_flags = titanlib::buddy_check(get_unflagged_points(), get_unflagged(values), get_unflagged(radius), get_unflagged(num_min), threshold, max_elev_diff, elev_gradient, min_std, num_iterations, get_unflagged(obs_to_check));
    merge(new_flags, indices);
}
void titanlib::Dataset::buddy_event_check(const vec& radius, const ivec& num_min, float event_threshold, float threshold, float max_elev_diff, float elev_gradient, int num_iterations, const ivec& obs_to_check, const ivec& indices) {
    ivec new_flags = titanlib::buddy_event_check(get_unflagged_points(), get_unflagged(values), get_unflagged(radius), get_unflagged(num_min), event_threshold, threshold, max_elev_diff, elev_gradient, num_iterations, get_unflagged(obs_to_check));
    merge(new_flags, indices);
}
void titanlib::Dataset::isolation_check(int num_min, float radius, float vertical_radius, const ivec& indices) {
    ivec new_flags = titanlib::isolation_check(get_unflagged_points(), num_min, radius, vertical_radius);
    merge(new_flags, indices);
}
void titanlib::Dataset::isolation_check(const ivec& num_min, const vec& radius, const vec& vertical_radius, const ivec& indices) {
    ivec new_flags = titanlib::isolation_check(get_unflagged_points(), get_unflagged(num_min), get_unflagged(radius), get_unflagged(vertical_radius));
    merge(new_flags, indices);
}
void titanlib::Dataset::duplicate_check(float radius, float vertical_range, const ivec& indices) {
    ivec new_flags = titanlib::duplicate_check(get_unflagged_points(), radius, vertical_range);
    merge(new_flags, indices);
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
    ivec new_flags = titanlib::metadata_check(get_unflagged_points(), check_lat, check_lon, check_elev, check_laf);
    merge_simple(new_flags, subset_valid(indices));
}

void titanlib::Dataset::merge(const ivec& new_flags, ivec indices) {
    if(indices.size() == 1 && indices[0] == -1) {
        indices.clear();
        indices.resize(flags.size());
        for(int i = 0; i < flags.size(); i++)
            indices[i] = i;
    }

    ivec flag_indices = get_unflagged_indices(); // These are the indices corresponding to each flag

    ivec indices0; // These are indices we want to update flags for
    indices0.reserve(indices.size());
    for(int i = 0; i < indices.size(); i++) {
        if(flags[indices[i]] == 0)
            indices0.push_back(indices[i]);
    }

    assert(new_flags.size() >= indices0.size());
    for(int i = 0; i < indices0.size(); i++) {
        int index = indices0[i];
        if(index >= flags.size() || index < 0)
            throw std::runtime_error("One or more indices are invalid");
        for(int j = 0; j < flag_indices.size(); j++) {
            if(flag_indices[j] == index)
                flags[index] = new_flags[j];
        }
    }
}

void titanlib::Dataset::merge_simple(const ivec& new_flags, ivec indices) {
    if(indices.size() == 1 && indices[0] == -1) {
        flags = new_flags;
        return;
    }

    std::cout << new_flags.size() << " " << indices.size() << std::endl;
    assert(new_flags.size() == indices.size());

    for(int i = 0; i < indices.size(); i++) {
        int index = indices[i];
        if(index >= flags.size() || index < 0)
            throw std::runtime_error("One or more indices are invalid");
        flags[index] = new_flags[i];
    }
}

ivec titanlib::Dataset::subset_valid(const ivec& indices) {
    if(indices.size() == 0)
        return ivec();
    ivec indices0 = indices;
    if(indices.size() == 1 && indices[0] == -1) {
        indices0.clear();
        indices0.resize(flags.size());
        for(int i = 0; i < flags.size(); i++)
            indices0[i] = i;
    }
    ivec output_indices;
    output_indices.reserve(indices0.size());
    for(int i = 0; i < indices0.size(); i++) {
        int index = indices0[i];
        if(flags[index] == 0) {
            output_indices.push_back(index);
        }
    }
    return output_indices;
}

ivec titanlib::Dataset::get_unflagged_indices() {
    ivec indices;
    indices.reserve(flags.size());
    for(int i = 0; i < flags.size(); i++) {
        if(flags[i] == 0) {
            indices.push_back(i);
        }
    }
    return indices;
}

Points titanlib::Dataset::get_unflagged_points() {
    ivec indices = get_unflagged_indices();
    vec ilats = points.get_lats();
    vec ilons = points.get_lons();
    vec ielevs = points.get_elevs();
    vec ilafs = points.get_lafs();
    vec lats(indices.size());
    vec lons(indices.size());
    vec elevs(indices.size());
    vec lafs(indices.size());
    for(int i = 0; i < indices.size(); i++) {
        int index = indices[i];
        assert(index < ilats.size());
        assert(index < ilons.size());
        assert(index < ielevs.size());
        assert(index < ilafs.size());
        lats[i] = ilats[index];
        lons[i] = ilons[index];
        elevs[i] = ielevs[index];
        lafs[i] = ilafs[index];
    }
    return Points(lats, lons, elevs, lafs);
}

Points titanlib::Dataset::subset_valid(const Points& input, const ivec& indices) {
    if(indices.size() == 0)
        return Points(vec(), vec(), vec(), vec());
    ivec indices0 = indices;
    if(indices.size() == 1 && indices[0] == -1) {
        indices0.clear();
        indices0.resize(flags.size());
        for(int i = 0; i < flags.size(); i++)
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
