#include <vector>
#include <algorithm>
#include <iostream>
#define _USE_MATH_DEFINES
#include <boost/accumulators/accumulators.hpp>
#include <boost/accumulators/statistics.hpp>
#include <boost/accumulators/statistics/variance.hpp>
#include <exception>
#include "titanlib.h"

using namespace titanlib;

ivec titanlib::buddy_event_check(const Points& points,
        const vec& values,
        const vec& radius,
        const ivec& num_min,
        float event_threshold,
        float threshold,
        float max_elev_diff,
        float elev_gradient,
        int num_iterations,
        const ivec& obs_to_check) {

    bool debug = false;
    const int s = values.size();
    // assert that the arrays we expect are of size s
    if(values.size() != s) {
        std::stringstream ss;
        ss << "Points (" << points.size() << ") and values (" << s << ") are not the same size";
        throw std::invalid_argument(ss.str());
    }
    if(radius.size() != s && radius.size() != 1) {
        std::stringstream ss;
        ss << "Radius (" << radius.size() << ") and values (" << s << ") are not the same size";
        throw std::invalid_argument(ss.str());
    }
    if(num_min.size() != s && num_min.size() != 1) {
        std::stringstream ss;
        ss << "'num_min' (" << num_min.size() << ") and values (" << s << ") are not the same size";
        throw std::invalid_argument(ss.str());
    }
    if(obs_to_check.size() != s && obs_to_check.size() != 1 && obs_to_check.size() !=0) {
        std::stringstream ss;
        ss << "'obs_to_check' (" << obs_to_check.size() << ") and values (" << s << ") are not the same size";
        throw std::invalid_argument(ss.str());
    }

    // Check that buddies min is more than 0
    for(int i = 0; i < num_min.size(); i++) {
        if(num_min[i] <= 0)
            throw std::runtime_error("Buddies_min must be > 0");
    }

    const vec& lats = points.get_lats();
    const vec& lons = points.get_lons();
    const vec& elevs = points.get_elevs();

    // resize the flags and set them to 0
    ivec flags(s, 0);
    for(int i = 0; i < num_min.size(); i++) {
        if(!titanlib::is_valid(values[i]))
            flags[i] = 1;
    }

    // if obs_to_check is empty then check all
    bool check_all = obs_to_check.size() != s;

    int num_removed_last_iteration = 0;
    for(int it = 0; it < num_iterations; it++) {
        ivec flags_prev = flags;
        #pragma omp parallel for
        for(int i = 0; i < values.size(); i++) {
            // is this one we are supposed to check?
            int b_i = (num_min.size() == s) ? i : 0;
            int d_i = (radius.size() == s) ? i : 0;
            if(flags_prev[i] != 0)
                continue;

            if(check_all || obs_to_check[i] == 1) {
                if(debug) {
                    std::cout << "point: " << lats[i] << " " << lons[i] << " " << elevs[i];
                    std::cout << ", and min buddies: " << num_min[b_i];
                    std::cout << '\n';
                }

                // get all neighbours that are close enough                 
                ivec neighbour_indices = points.get_neighbours(lats[i], lons[i], radius[d_i], false);

                int n_buddies = 0;
                std::vector<bool> buddy_events;
                // based on tree do have enough neighbours? 
                if(neighbour_indices.size() >= num_min[b_i]) {
                    // loop over everything that was near enough
                    // count buddies and make list of values (adjusting for height diff if needed)
                    for(int j = 0; j < neighbour_indices.size(); j++) {
                        if(flags_prev[neighbour_indices[j]] != 0)
                            continue;

                        // don't use ones that differ too much in height (max_elev_diff)
                        if(max_elev_diff > 0) {
                            float elev_diff = fabs(elevs[neighbour_indices[j]] - elevs[i]);
                            if(elev_diff <= max_elev_diff) {
                                // correction for the elevation differences (add or subtract -0.0065 degC/m)
                                // m difference from point in question
                                float elev_diff = elevs[i] - elevs[neighbour_indices[j]];
                                //std::cout << "height diff: " << elev_diff;
                                float adjusted_value = values[neighbour_indices[j]] + (elev_diff * elev_gradient);
                                //std::cout << ", adjusted value: " << adjusted_value;
                                //std::cout << '\n';
                                bool is_event = adjusted_value < event_threshold;

                                buddy_events.push_back(is_event);
                                n_buddies++;
                            }
                            else {
                                if(debug) {
                                    std::cout << "too much height difference: " << elev_diff << '\n';
                                }
                            }
                        }
                        // if max_elev_diff is negative then don't check elevation difference
                        else {
                            // can use this station
                            int is_event = values[neighbour_indices[j]] < event_threshold;
                            buddy_events.push_back(is_event);
                            n_buddies++;
                        }
                    }

                }
                if(debug) {
                    std::cout << "buddies: " << n_buddies << '\n';
                }
                if(n_buddies >= num_min[b_i]) {
                    // compute the average and standard deviation of the values
                    int count_event = 0;
                    for(int k = 0; k < buddy_events.size(); k++) {
                        if(buddy_events[k] == 1)
                            count_event++;
                    }
                    float fraction_event = float(count_event) / buddy_events.size();
                    bool is_event = values[i] < event_threshold;

                    if(threshold < 1) {
                        if(is_event && fraction_event <= threshold)
                            flags[i] = 1;
                        else if(!is_event && (1 - fraction_event) <= threshold)
                            flags[i] = 1;
                    }
                    else {
                        if(is_event && count_event <= threshold)
                            flags[i] = 1;
                        else if(!is_event && (buddy_events.size() - count_event <= threshold))
                            flags[i] = 1;
                    }
                    if(debug && flags[i] == 1) {
                        std::cout << "value: " << values[i] << '\n';
                        std::cout << "is_event: " << is_event << '\n';
                        std::cout << "event_threshold: " << event_threshold << '\n';
                        std::cout << "threshold: " << threshold << '\n';
                        std::cout << "count: " << buddy_events.size() << '\n';
                        std::cout << "count_event: " << count_event << '\n';
                        std::cout << "fraction_event: " << fraction_event << '\n';
                        std::cout << "flag: " << flags[i] << '\n';
                        std::cout << '\n';
                    }
                }
            }
        }
        // Check if we need to stop early
        int num_removed = 0;
        for(int i = 0; i < s; i++) {
            if(flags[i] != 0)
                num_removed++;
        }
        int num_removed_curr_iteration = num_removed - num_removed_last_iteration;
        if(debug) {
            std::cout << "iteration, number of bad observations: " << it + 1 << ", " << num_removed_curr_iteration << '\n';
        }
        if(num_removed_curr_iteration == 0) {
            if(debug)
                std::cout << "Stopping early after iteration " << it + 1 << std::endl;
            break;
        }
        num_removed_last_iteration = num_removed_curr_iteration;
    }

    return flags;
}
