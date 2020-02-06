#ifndef TITANLIB_H
#define TITANLIB_H
#include <iostream>
#include <vector>
typedef std::vector<float> fvec;
typedef std::vector<int> ivec;

/** Titanlib
*/
namespace titanlib {
    float calc_gamma(float shape, float scale);

    /** Spatial Consistency Test
      * @param lats vector of latitudes
      * @param flags output vector of flags
      */
    bool sct(const fvec lats,
            const fvec lons,
            const fvec elevs,
            const fvec values,
            int nmin,
            int nmax,
            int nminprof,
            float dzmin,
            float dhmin,
            float dz,
            const fvec t2pos,
            const fvec t2neg,
            const fvec eps2,
            ivec& flags);


    bool range_check(const fvec values,
            const fvec min,
            const fvec max,
            ivec& flags);

    bool range_check_climatology(const fvec lats,
            const fvec lons,
            const fvec elevs,
            const fvec values,
            const fvec min,
            const fvec max,
            ivec& flags);

}
#endif
