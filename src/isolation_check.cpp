#include <vector>
#include <math.h>
#include "titanlib.h"
#include <assert.h>
#include <iostream>

ivec titanlib::isolation_check(const fvec& lats,
        const fvec& lons,
        int num_min,
        float radius) {
    titanlib::KDTree tree(lats, lons);
    ivec flags(lats.size(), 0);

    for(int i = 0; i < lats.size(); i++) {
        int num = tree.get_num_neighbours(lats[i], lons[i], radius, 0, false);
        // std::cout << i << " " << num << std::endl;
        if(num < num_min) {
            flags[i] = 1;
        }
    }

    return flags;

}

ivec titanlib::isolation_check(const fvec& lats,
        const fvec& lons,
        const fvec& elevs,
        int num_min,
        float radius,
        float dz) {

    titanlib::KDTree tree(lats, lons);
    ivec flags(lats.size(), 0);

    for(int i = 0; i < lats.size(); i++) {
        ivec indices = tree.get_neighbours(lats[i], lons[i], radius, 0, false);
        int num = 0;
        for(int j = 0; j < indices.size(); j++) {
            int index = indices[j];
            if(fabs(elevs[index] - elevs[i]) < dz)
                num++;
        }
        // std::cout << i << " " << num << std::endl;
        if(num < num_min ) {
            flags[i] = 1;
        }
    }

    return flags;
}
