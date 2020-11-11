#include "titanlib.h"
#include <vector>
#include <iostream>

using namespace titanlib;

ivec titanlib::duplicate_check(const Points& points, float radius) {

    const vec& lats = points.get_lats();
    const vec& lons = points.get_lons();
    const vec& elevs = points.get_elevs();

    ivec keep;
    ivec checked(points.size(), 0); // Which points have been checked
    keep.reserve(points.size());

    for(int i = 0; i < points.size(); i++) {
        if(checked[i] == 1)
            continue;
        ivec indices = points.get_neighbours(lats[i], lons[i], radius, true);
        int num = indices.size();
        checked[i] = 1;
        keep.push_back(i);
        if(num > 1) {
            for(int j = 0; j < num; j++) {
                checked[indices[j]] = 1;
            }
        }
    }
    ivec flags(points.size(), 1);
    for(int i = 0; i < keep.size(); i++) {
        flags[keep[i]] = 0;
    }

    return flags;
}
