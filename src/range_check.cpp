#include <vector>
#include "titanlib.h"
#include <assert.h>


bool titanlib::range_check(const fvec lats,
        const fvec lons,
        const fvec elevs,
        const fvec values,
        const fvec min,
        const fvec max,
        ivec& flags) {

    // loop over all the lats/lons/elevs + value 
    // either min/max has length 1 or is the same length as the other vecs
    const int s = lats.size();
    if( lons.size() != s || elevs.size() != s || values.size() != s ) { return false; }
    if( (min.size() != s && min.size() != 1) || (max.size() != s && max.size() != 1) ) { return false; }

    flags.resize(s, 0);

    for(int i = 0; i < s; i++) {
        // leave the index to 0 if its the same max/min applied to everything
        // else same as loop
        int min_i = (min.size() == s) ? i : 0;
        int max_i = (max.size() == s) ? i : 0;

        // loop over the vectors and set the flags (0 = ok and 1 = bad)
        if(values[i] < min[min_i] || values[i] > max[max_i]) {
            flags[i] = 1;
        }
    }

    return true;

}
