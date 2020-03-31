#include <vector>
#include <algorithm>
#include <iostream>
#define _USE_MATH_DEFINES
#include <boost/accumulators/accumulators.hpp>
#include <boost/accumulators/statistics.hpp>
#include <boost/accumulators/statistics/variance.hpp>
#include "titanlib.h"

int titanlib::buddy_check(const fvec& lats,
        const fvec& lons,
        const fvec& elevs,
        const fvec& values,
        const fvec& radius,
        const ivec& buddies_min,
        const fvec& thresholds,
        float diff_elev_max,
        float elev_gradient,
        float min_std,
        int num_iterations,
        ivec& flags,
        const ivec obs_to_check) {

    bool debug = false;
    const int s = values.size();
    // assert that the arrays we expect are of size s
    if( lats.size() != s || lons.size() != s || elevs.size() != s || values.size() != s) { return 1; }
    if( radius.size() != s && radius.size() != 1 ) { return 1; }
    if( (buddies_min.size() != s && buddies_min.size() != 1) || (thresholds.size() != s && thresholds.size() != 1) ) { return 1; }
    if( (obs_to_check.size() != s && obs_to_check.size() != 1 && obs_to_check.size() !=0) ) { return 1; }

    // Check that buddies min is more than 0
    for(int i = 0; i < buddies_min.size(); i++) {
        if(buddies_min[i] <= 0)
            return 1;
    }

    // create the KD tree to be used later
    titanlib::KDTree tree(lats, lons);
    // resize the flags and set them to 0
    flags.resize(s, 0);
    // if obs_to_check is empty then check all
    bool check_all = (obs_to_check.size() == s) ? false : true;

    // loop over all the observations
    for(int it = 0; it < num_iterations; it++) {
        for(int i = 0; i < values.size(); i++) {
            // is this one we are supposed to check?
            int b_i = (buddies_min.size() == s) ? i : 0;
            int d_i = (radius.size() == s) ? i : 0;
            int t_i = (thresholds.size() == s) ? i : 0;
            if(flags[i] != 0)
                continue;
            if( ((!check_all && obs_to_check[i] == 1) || check_all) ) {
                if(debug) {
                    std::cout << "point: " << lats[i] << " " << lons[i] << " " << elevs[i];
                    std::cout << ", and min buddies: " << buddies_min[b_i];
                    std::cout << '\n';
                }

                // get all neighbours that are close enough
                ivec neighbour_indices = tree.get_neighbours(lats[i], lons[i], radius[d_i], 0, false);

                int n_buddies = 0;
                fvec list_buddies;
                // based on tree do have enough neighbours? 
                if(neighbour_indices.size() > buddies_min[b_i]) {
                    // loop over everything that was near enough
                    // count buddies and make list of values (adjusting for height diff if needed)
                    for(int j = 0; j < neighbour_indices.size(); j++) {
                        // don't use ones that differ too much in height (diff_elev_max)
                        if(flags[neighbour_indices[j]] == 0) {
                            if(diff_elev_max > 0) {
                                float elev_diff = fabs(elevs[neighbour_indices[j]] - elevs[i]);
                                if(elev_diff <= diff_elev_max) {
                                    // correction for the elevation differences (add or subtract -0.0065 degC/m)
                                    // m difference from point in question
                                    float elev_diff = elevs[neighbour_indices[j]] - elevs[i];
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
                            // if diff_elev_max is negative then don't check elevation difference
                            else {
                                // can use this station
                                list_buddies.push_back(values[neighbour_indices[j]]);
                                n_buddies++;
                            }
                        }
                    }

                }
                if(debug) {
                    std::cout << "buddies: " << n_buddies << '\n';
                }
                if(n_buddies >= buddies_min[b_i]) {
                    // compute the average and standard deviation of the values
                    boost::accumulators::accumulator_set<float, boost::accumulators::features<boost::accumulators::tag::mean, boost::accumulators::tag::variance>> acc;
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
                    if(std < min_std) {
                        std = min_std;
                    }
                    float pog = fabs(values[i] - mean)/std;
                    if(pog > thresholds[t_i]) {
                        flags[i] = 1;
                    }
                }
            }
        }
    }

    return 0;
}
