#include "titanlib.h"
#include <vector>
#include <iostream>

using namespace titanlib;

ivec titanlib::duplicate_check(const Points& points, float radius, float vertical_range) {

    const vec& lats = points.get_lats();
    const vec& lons = points.get_lons();
    const vec& elevs = points.get_elevs();

    ivec keep;
    ivec checked(points.size(), 0); // Which points have been checked
    keep.reserve(points.size());
    bool do_elevation_check = titanlib::is_valid(vertical_range);
    if(do_elevation_check) {
        for(int i = 0; i < points.size(); i++) {
            if(!titanlib::is_valid(elevs[i]))
                checked[i] = 1;
        }
    }

    // TODO: Can this be parallelized? Looks difficult...
    for(int i = 0; i < points.size(); i++) {
        if(checked[i] == 1)
            continue;
        ivec indices = points.get_neighbours(lats[i], lons[i], radius, true);
        const int num = indices.size();
        checked[i] = 1;
        keep.push_back(i);
        if(num > 1) {
            // Count how many are within the elevation range
            if(do_elevation_check) {
                for(int j = 0; j < indices.size(); j++) {
                    float curr_elev = elevs[indices[j]];
                    if(!titanlib::is_valid(curr_elev))
                        checked[indices[j]] = 1;
                    else if(titanlib::is_valid(curr_elev) && (fabs(elevs[i] - curr_elev) <= vertical_range))
                        checked[indices[j]] = 1;
                }
            }
            else {
                for(int j = 0; j < num; j++) {
                    checked[indices[j]] = 1;
                }
            }
        }
    }
    ivec flags(points.size(), 1);
    for(int i = 0; i < keep.size(); i++) {
        flags[keep[i]] = 0;
    }

    return flags;
}
