#include <vector>
#include <algorithm>
#include <iostream>
#define _USE_MATH_DEFINES
#include <math.h>
#include <boost/accumulators/accumulators.hpp>
#include <boost/accumulators/statistics.hpp>
#include <boost/accumulators/statistics/variance.hpp>
#include "titanlib.h"

float calc_distance(float lat1, float lon1, float lat2, float lon2);
float deg2rad(float deg);

bool titanlib::buddy_check(const fvec lats,
        const fvec lons,
        const fvec elevs,
        const fvec values,
        const ivec obs_to_check, // 1 = true, 0 == false 
        const fvec distance_lim,
        const fvec priorities, // need to know the range of priorities?
        const fvec buddies_min,
        const fvec thresholds,
        float diff_elev_max,
        bool adjust_for_elev_diff,
        ivec& flags) {

    const int s = values.size();
    // assert that the arrays we expect are of size s
    if( lats.size() != s || lons.size() != s || elevs.size() != s || priorities.size() != s || obs_to_check.size() != s) { return false; }
    if( distance_lim.size() != s || priorities.size() != s || buddies_min.size() != s || thresholds.size() != s ) { return false; }
    // resize the flags and set them to 0
    flags.resize(s, 0);

    // loop to find all levels of priorities
    fvec p_list; 
    for(int x = 0; x < priorities.size(); x++) {
        std::cout << priorities[x] << " ";
    }
    std::cout << '\n';
    for(int p = 0; p < priorities.size(); p++) {
        std::cout << "priority: " << priorities[p];
        std::cout << '\n';
        std::vector<float>::iterator it;
        it = std::find(p_list.begin(),p_list.end(),priorities[p]);
        if(it != p_list.end()) {
            std::cout << "element already in priorities list: " << priorities[p];
            std::cout << '\n';
        }
        else {
            std::cout << "adding element to priorities list: " << priorities[p];
            std::cout << '\n';
            p_list.push_back(priorities[p]);
        }
    }
    for(int x = 0; x < priorities.size(); x++) {
        std::cout << priorities[x] << " ";
    }
    std::cout << '\n';
    std::sort(p_list.begin(), p_list.end());
    std::cout << "priorities list sorted: ";
    for(int x = 0; x < p_list.size(); x++) {
        std::cout << p_list[x] << " ";
    }
    std::cout << '\n';

    // for each priority
    for(int p = 0; p < p_list.size(); p++) {
        // loop over all the observations
        for(int i = 0; i < values.size(); i++) {
            // is this one we are supposed to check?
            // and is it of the priority we are currently looking at?
            if(obs_to_check[i] == 1 && priorities[i] == p_list[p]) {
                std::cout << "point: " << lats[i] << " " << lons[i] << " " << elevs[i]; 
                std::cout << ", priority: " << priorities[i];
                std::cout << ", and min buddies: " << buddies_min[i];
                std::cout << '\n';

                // loop through all the obs and check if they are close enough
                // keep track of how many buddies it has, to ensure it is more than buddies_min 
                int n_buddies = 0;
                fvec list_buddies;
                for(int j = 0; j < values.size(); j++) { 
                    // do not compare it to itself
                    if(j != i) {
                        float d = calc_distance(lats[i], lons[i], lats[j], lons[j]);
                        std::cout << "other point: " << lats[j] << " " << lons[j] << " " << elevs[j];
                        std::cout << ", distance from 1st point: " << d;
                        std::cout << '\n';
                        // find all the observations within the prescribed distance
                        if(d <= distance_lim[i]) {
                            // don't use ones that differ too much in height (diff_elev_max)
                            if(diff_elev_max > 0) {
                                float elev_diff = fabs(elevs[j] - elevs[i]);
                                //std::cout << "elev diff: " << elev_diff << '\n';
                                if(elev_diff <= diff_elev_max) {
                                    // can use this station
                                    if(adjust_for_elev_diff) {
                                        // correction for the elevation differences (add or subtract -0.0065 degC/m)
                                        // m difference from point in question
                                        float elev_diff = elevs[j] - elevs[i];
                                        std::cout << "height diff: " << elev_diff;
                                        float adjusted_value = values[j] + (elev_diff * 0.0065);
                                        std::cout << ", adjusted value: " << adjusted_value;
                                        std::cout << '\n';
                                        list_buddies.push_back(adjusted_value);
                                    }
                                    else {
                                        list_buddies.push_back(values[j]);
                                    }
                                    n_buddies++;
                                }
                                else {
                                    std::cout << "too much height difference!" << '\n';
                                }
                            }
                            // if diff_elev_max is negative then don't check elevation difference
                            else {
                                // can use this station
                                list_buddies.push_back(values[j]);
                                n_buddies++;
                            }
                        }
                        else {
                            std::cout << "too far between points!" << '\n';
                        }
                    }
                }
                std::cout << "buddies: " << n_buddies << '\n';
                if(n_buddies >= buddies_min[i]) {
                    // compute the average and standard deviation of the values
                    boost::accumulators::accumulator_set<float, boost::accumulators::features<boost::accumulators::tag::mean, boost::accumulators::tag::variance>> acc;
                    for(int k = 0; k < list_buddies.size(); k++) {
                        acc(list_buddies[k]);
                    }
                    float mean = boost::accumulators::mean(acc);
                    float variance = boost::accumulators::variance(acc);
                    std::cout << "mean: " << mean << '\n';
                    std::cout << "variance: " << variance << '\n';

                    float pog = fabs(values[i] - mean);
                    if(pog > thresholds[i]) {
                        flags[i] = 1;
                    }
                }
            }
        }
    }

    return true;
}

float calc_distance(float lat1, float lon1, float lat2, float lon2) {
    if(!(fabs(lat1) <= 90 && fabs(lat2) <= 90 && fabs(lon1) <= 360 && fabs(lon2) <= 360)) {
        std::cout << " Cannot calculate distance, invalid lat/lon: (" << lat1 << "," << lon1 << ") (" << lat2 << "," << lon2 << ")";
        std::cout << '\n';
    }
    if(lat1 == lat2 && lon1 == lon2)
        return 0;

    double lat1r = deg2rad(lat1);
    double lat2r = deg2rad(lat2);
    double lon1r = deg2rad(lon1);
    double lon2r = deg2rad(lon2);
    double radiusEarth = 6.378137e6;

    double ratio = cos(lat1r)*cos(lon1r)*cos(lat2r)*cos(lon2r)
                   + cos(lat1r)*sin(lon1r)*cos(lat2r)*sin(lon2r)
                   + sin(lat1r)*sin(lat2r);
    double dist = acos(ratio)*radiusEarth;
    return (float) dist;
}

float deg2rad(float deg) {
   return (deg * M_PI / 180);
}

