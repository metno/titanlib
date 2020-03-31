#include <vector>
#include <math.h>
#include "titanlib.h"
#include <assert.h>
#include <iostream>
extern "C" {
#include "sct_smart_boxes.h"
}

int titanlib::sct(const fvec& lats,
        const fvec& lons,
        const fvec& elevs,
        const fvec& values,
        int nmin,
        int nmax,
        int nminprof,
        float dzmin,
        float dhmin,
        float dz,
        const fvec& t2pos,
        const fvec& t2neg,
        const fvec& eps2,
        fvec& sct,
        fvec& rep,
        ivec& flags) {

    // Check inputs
    if(nmin == 0 || nmax == 0)
        return 1;
    fvec x, y;
    titanlib::util::convert_to_proj(lats, lons, "+proj=lcc +lat_0=63 +lon_0=15 +lat_1=63 +lat_2=63 +no_defs +R=6.371e+06", x, y);

    dvec dx(x.begin(), x.end());
    dvec dy(y.begin(), y.end());

    dvec delevs(elevs.begin(), elevs.end());
    dvec dvalues(values.begin(), values.end());
    dvec dt2pos(t2pos.begin(), t2pos.end());
    dvec dt2neg(t2neg.begin(), t2neg.end());
    dvec deps2(eps2.begin(), eps2.end());
    int N = lats.size();

    double ddzmin = dzmin;
    double ddhmin = dhmin;
    double ddz = dz;
    // int nmax = 200;
    // int nmin = 50;
    // int nminprof = 10;
    // double dzmin = 100;
    // double dhmin = 10000;
    // double dz = 30;
    dvec dsct;
    dvec drep;
    ivec boxids;
    boxids.resize(N, 0);
    flags.resize(N, 0);
    dsct.resize(N, 0);
    drep.resize(N, 0);

    sct_smart_boxes(&N, &dx[0], &dy[0], &delevs[0], &dvalues[0], &nmax, &nmin, &nminprof, &ddzmin, &ddhmin, &ddz, &dt2pos[0], &dt2neg[0], &deps2[0], &flags[0], &dsct[0], &drep[0], &boxids[0]);
    sct = fvec(dsct.begin(), dsct.end());
    rep = fvec(drep.begin(), drep.end());

    return 0;
}
