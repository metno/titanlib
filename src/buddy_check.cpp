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

ivec titanlib::buddy_check(const Points& points,
        const vec& values,
        const vec& radius,
        const ivec& num_min,
        float threshold,
        float max_elev_diff,
        float elev_gradient,
        float min_std,
        int num_iterations,
        const ivec& obs_to_check) {

    bool debug = false;
    const int s = points.size();
    // assert that the arrays we expect are of size s
    if(values.size() != s) {
        std::stringstream ss;
        ss <<"'values' (" << values.size() << ") does not have the same length as points (" << s << ")" << std::endl;
        throw std::invalid_argument(ss.str());
    }
    else if(radius.size() != s && radius.size() != 1) {
        std::stringstream ss;
        ss <<"'radius' (" << radius.size() << ") does not have the same length as points (" << s << ")" << std::endl;
        throw std::invalid_argument(ss.str());
    }
    if(num_min.size() != s && num_min.size() != 1) {
        std::stringstream ss;
        ss <<"'num_min' (" << num_min.size() << ") does not have the same length as points (" << s << ")" << std::endl;
        throw std::invalid_argument(ss.str());
    }
    if(obs_to_check.size() != s && obs_to_check.size() != 1 && obs_to_check.size() !=0) {
        std::stringstream ss;
        ss <<"'obs_to_check' (" << obs_to_check.size() << ") does not have the same length as points (" << s << ")" << std::endl;
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
    if(s==0) {
        // if no data return empty array
        return flags;
    }
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
        for(int i = 0; i < s; i++) {
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
                vec list_buddies;
                // based on tree do have enough neighbours? 
                if(neighbour_indices.size() >= num_min[b_i]) {
                    // loop over everything that was near enough
                    // count buddies and make list of values (adjusting for height diff if needed)
                    for(int j = 0; j < neighbour_indices.size(); j++) {
                        // don't use ones that differ too much in height (max_elev_diff)
                        if(flags_prev[neighbour_indices[j]] != 0)
                            continue;

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
                                list_buddies.push_back(adjusted_value);
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
                            list_buddies.push_back(values[neighbour_indices[j]]);
                            n_buddies++;
                        }
                    }

                }
                if(debug) {
                    std::cout << "buddies: " << n_buddies << '\n';
                }
                if(n_buddies >= num_min[b_i]) {
                    // compute the average and standard deviation of the values
                    boost::accumulators::accumulator_set<float, boost::accumulators::features<boost::accumulators::tag::mean, boost::accumulators::tag::variance> > acc;
                    for(int k = 0; k < list_buddies.size(); k++) {
                        acc(list_buddies[k]);
                    }
                    float mean = boost::accumulators::mean(acc);
                    float variance = boost::accumulators::variance(acc);
                    if(debug) {
                        std::cout << "mean: " << mean << '\n';
                        std::cout << "variance: " << variance << '\n';
                    }

                    float std = sqrt(variance);
                    float std_adjusted = sqrt(variance + variance / n_buddies);
                    if(std_adjusted < min_std) {
                        std_adjusted = min_std;
                    }
                    float pog = fabs(values[i] - mean)/std_adjusted;
                    if(pog > threshold) {
                        flags[i] = 1;
                    }
                }
            } // end if observation is to check
        } // end loop over values

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
    } // end loop over num_iterations

    return flags;
}
