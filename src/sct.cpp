#include <vector>
#include <math.h>
#include "titanlib.h"
#include <assert.h>
#include <iostream>
extern "C" {
#include "sct_smart_boxes.h"
}

bool titanlib::sct(const fvec lats,
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
        ivec& flags) {

    dvec dlats(lats.begin(), lats.end());
    dvec dlons(lons.begin(), lons.end());
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
    dvec sct;
    dvec rep;
    ivec boxids;
    sct.resize(N, 0);
    rep.resize(N, 0);
    boxids.resize(N, 0);
    flags.resize(N, 0);

    sct_smart_boxes(&N, &dlons[0], &dlats[0], &delevs[0], &dvalues[0], &nmax, &nmin, &nminprof, &ddzmin, &ddhmin, &ddz, &dt2pos[0], &dt2neg[0], &deps2[0], &flags[0], &sct[0], &rep[0], &boxids[0]);

    return true;
}
