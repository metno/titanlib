#include <vector>
#include <math.h>
#include "titanlib.h"
#include <assert.h>
#include <iostream>

using namespace titanlib;

ivec titanlib::isolation_check(const vec& lats,
        const vec& lons,
        int num_min,
        float radius) {
    titanlib::KDTree tree(lats, lons);
    ivec flags(lats.size(), 0);

    for(int i = 0; i < lats.size(); i++) {
        int num = tree.get_num_neighbours(lats[i], lons[i], radius, false);
        // std::cout << i << " " << num << std::endl;
        if(num < num_min) {
            flags[i] = 1;
        }
    }

    return flags;

}

ivec titanlib::isolation_check(const vec& lats,
        const vec& lons,
        const vec& elevs,
        int num_min,
        float radius,
        float vertical_radius) {

    titanlib::KDTree tree(lats, lons);
    ivec flags(lats.size(), 0);

    for(int i = 0; i < lats.size(); i++) {
        ivec indices = tree.get_neighbours(lats[i], lons[i], radius, false);
        int num = 0;
        for(int j = 0; j < indices.size(); j++) {
            int index = indices[j];
            if(fabs(elevs[index] - elevs[i]) < vertical_radius)
                num++;
        }
        // std::cout << i << " " << num << std::endl;
        if(num < num_min ) {
            flags[i] = 1;
        }
    }

    return flags;
}
