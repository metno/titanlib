#include <vector>
#include <math.h>
#include "titanlib.h"
#include <assert.h>
#include <iostream>

using namespace titanlib;

ivec titanlib::isolation_check(const Points& points,
        int num_min,
        float radius,
        float vertical_radius) {

    const int s = points.size();
    ivec flags(s, 0);
    const vec& lats = points.get_lats();
    const vec& lons = points.get_lons();
    const vec& elevs = points.get_elevs();

    ivec num_min_array(s);
    vec radius_array(s);
    vec vertical_radius_array(s);
    for(int i = 0; i < s; i++) {
        num_min_array[i] = num_min;
        radius_array[i] = radius;
        vertical_radius_array[i] = vertical_radius;
    }

    return isolation_check(points, num_min_array, radius_array, vertical_radius_array);
}

ivec titanlib::isolation_check(const Points& points,
        const ivec& num_min,
        const vec& radius,
        const vec& vertical_radius) {

    const int s = points.size();
    if(num_min.size() != s) {
        std::stringstream ss;
        ss <<"'num_min' (" << num_min.size() << ") does not have the same length as points (" << s << ")" << std::endl;
        throw std::invalid_argument(ss.str());
    }
    else if(radius.size() != s) {
        std::stringstream ss;
        ss <<"'radius' (" << radius.size() << ") does not have the same length as points (" << s << ")" << std::endl;
        throw std::invalid_argument(ss.str());
    }
    else if(vertical_radius.size() > 0 && vertical_radius.size() != s) {
        std::stringstream ss;
        ss <<"'vertical_radius' (" << vertical_radius.size() << ") does not have the same length as points (" << s << ")" << std::endl;
        throw std::invalid_argument(ss.str());
    }
    ivec flags(s, 0);
    const vec& lats = points.get_lats();
    const vec& lons = points.get_lons();
    const vec& elevs = points.get_elevs();

    #pragma omp parallel for
    for(int i = 0; i < s; i++) {
        if(!titanlib::is_valid(lats[i]) || !titanlib::is_valid(lons[i])) {
            flags[i] = 1;
            continue;
        }
        if(vertical_radius.size() > 0 && titanlib::is_valid(vertical_radius[i])) {
            if(!titanlib::is_valid(elevs[i])) {
                flags[i] = 1;
                continue;
            }
            ivec indices = points.get_neighbours(lats[i], lons[i], radius[i], false);
            int num = 0;
            for(int j = 0; j < indices.size(); j++) {
                int index = indices[j];
                if(titanlib::is_valid(elevs[index]) && titanlib::is_valid(elevs[i])) {
                    if(fabs(elevs[index] - elevs[i]) <= vertical_radius[i])
                        num++;
                }
            }
            if(num < num_min[i]) {
                flags[i] = 1;
            }
        }
        else {
            // Faster version when we don't need to check elevations
            int num = points.get_num_neighbours(lats[i], lons[i], radius[i], false);
            if(num < num_min[i]) {
                flags[i] = 1;
            }
        }
    }

    return flags;
}
