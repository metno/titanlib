#include <vector>
#include <math.h>
#include "titanlib.h"
#include <assert.h>
#include <iostream>
extern "C" {
#include "sct_smart_boxes.h"
}

int titanlib::sct_window(const fvec& lats,
            const fvec& lons,
            const fvec& elevs,
            const fvec& values,
            int nminprof,
            float radius,
            float dzmin,
            float dhmin,
            float dz,
            const fvec& pos,
            const fvec& neg,
            const fvec& eps2,
            fvec& sct,
            ivec& flags) {

    const int s = values.size();
    // assert that the arrays we expect are of size s
    if( lats.size() != s || lons.size() != s || elevs.size() != s || values.size() != s) { return 1; }

    // create the KD tree to be used later
    titanlib::KDTree tree(lats, lons);

    // resize the flags and set them to 0
    flags.resize(s, 0);

    for(int i = 0; i < values.size(); i++) {
        // get all neighbours that are close enough to this point
        ivec neighbour_indices = tree.get_neighbours(lats[i], lons[i], radius, 0, false);
        int s_box = neighbour_indices.size()+1;
        fvec x(s_box);
        fvec y(s_box);
        fvec z(s_box);
        fvec e(s_box);
        fvec t(s_box);
        // add the observation itself 
        titanlib::util::convert_coordinates(lats[i], lons[i], x[0], y[0], z[0]);
        e[0] = elevs[i];
        t[0] = values[i];
        // add all the obs that are close enough 
        for(int j = 1; s_box; j++) {
            titanlib::util::convert_coordinates(lats[neighbour_indices[j]], lons[neighbour_indices[j]], x[j], y[j], z[j]);
            e[j] = elevs[neighbour_indices[j]];
            t[j] = values[neighbour_indices[j]];
        }
        ivec obs_check;
        obs_check.resize(s_box, 0);
        // just check the first one
        obs_check[0] = 1;

        // fvec x2, y2;
        // titanlib::util::convert_to_proj(x, y, "+proj=lcc +lat_0=63 +lon_0=15 +lat_1=63 +lat_2=63 +no_defs +R=6.371e+06", x, y);

        dvec d_x(x.begin(), x.end());
        dvec d_y(y.begin(), y.end());
        dvec d_z(z.begin(), z.end());

        dvec d_elevs(e.begin(), e.end());
        dvec d_values(t.begin(), t.end());
        dvec d_t2pos(pos.begin(), pos.end());
        dvec d_t2neg(neg.begin(), neg.end());
        dvec d_eps2(eps2.begin(), eps2.end());
        int N = x.size();

        ivec d_check(obs_check.begin(), obs_check.end());

        double d_dzmin = dzmin;
        double d_dhmin = dhmin;
        double d_dz = dz;

        dvec d_sct;
        dvec rep;
        ivec boxids;
        ivec flags_local;
        sct.resize(N, 0);
        rep.resize(N, 0);
        boxids.resize(N, 0);
        flags_local.resize(N, 0);
        d_sct.resize(N, 0);

        spatial_consistency_test_mod(&N, &d_check[0], &d_x[0], &d_y[0], &d_z[0], &d_elevs[0], &d_values[0], &nminprof, &d_dzmin, &d_dhmin, &d_dz, &d_t2pos[0], &d_t2neg[0], &d_eps2[0], &flags_local[0], &d_sct[0], &rep[0]);
        sct = fvec(d_sct.begin(), d_sct.end());

        // did it flag any of the obs we wanted to check?
        for(int k = 0; k < obs_check.size(); k++) {
            if(obs_check[k] == 0) {
                if(k == 0)
                    flags[i] = flags_local[k];
                else
                    flags[neighbour_indices[k]] = flags_local[k];
            }
        }
    }

    return 0;
}
