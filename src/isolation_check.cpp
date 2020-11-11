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

    ivec flags(points.size(), 0);
    const vec& lats = points.get_lats();
    const vec& lons = points.get_lons();
    const vec& elevs = points.get_elevs();

    for(int i = 0; i < points.size(); i++) {
        if(titanlib::is_valid(vertical_radius)) {
            ivec indices = points.get_neighbours(lats[i], lons[i], radius, false);
            int num = 0;
            for(int j = 0; j < indices.size(); j++) {
                int index = indices[j];
                if(fabs(elevs[index] - elevs[i]) < vertical_radius)
                    num++;
            }
            if(num < num_min ) {
                flags[i] = 1;
            }
        }
        else {
            // Faster version when we don't need to check elevations
            int num = points.get_num_neighbours(lats[i], lons[i], radius, false);
            if(num < num_min) {
                flags[i] = 1;
            }
        }
    }

    return flags;
}
