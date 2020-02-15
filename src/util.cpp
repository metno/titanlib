#include <cmath>
#include "titanlib.h"
#include <libalglib/interpolation.h>

bool titanlib::util::convert_coordinates(const fvec& lats, const fvec& lons, fvec& x_coords, fvec& y_coords, fvec& z_coords) {
    int N = lats.size();
    x_coords.resize(N);
    y_coords.resize(N);
    z_coords.resize(N);
    for(int i = 0; i < N; i++) {
        convert_coordinates(lats[i], lons[i], x_coords[i], y_coords[i], z_coords[i]);
    }

    return true;
}

bool titanlib::util::convert_coordinates(float lat, float lon, float& x_coord, float& y_coord, float& z_coord) {

    float earth_radius = 6.37e6;
    double lonr = M_PI / 180 * lon;
    double latr = M_PI / 180 * lat;
    // std::cout << lon << " " << lat << std::endl;
    x_coord = std::cos(latr) * std::cos(lonr) * earth_radius;
    y_coord = std::cos(latr) * std::sin(lonr) * earth_radius;
    z_coord = std::sin(latr) * earth_radius;
    return true;
}
/*
ivec titanlib::util::nearest_neighbours(const fvec& lats, const fvec& lons, float radius, float lat, float lon) {
    fvec x, y, z;
    int nx = 3;
    int ny = 1;
    int normtype = 2;
    convert_coordinates(lats, lons, x, y, z);
    assert(x.size() == y.size() == z.size());
    alglib::real_2d_array a; // = "[[0,0,0],[0,1,0],[1,0,0],[1,1,0]]";
    a.setlength(lats.size(), 4);
    for(int i = 0; i < lats.size(); i++) {
        a[i][0] = x[i];
        a[i][1] = y[i];
        a[i][2] = z[i];
        a[i][3] = i;
    }
    alglib::kdtree tree;

    alglib::real_1d_array b;
    fvec x0, y0, z0;
    fvec lat0, lon0;
    lat0.resize(1, lat);
    lon0.resize(1, lon);
    titanlib::convert_coordinates(lat0, lon0, x0, y0, z0);
    b.setlength(3);
    b[0] = x0[0];
    b[1] = y0[0];
    b[2] = z0[0];

    alglib::kdtreebuild(a, nx, ny, normtype, tree);
    std::cout << 1 << std::endl;
    int num = alglib::kdtreequeryrnn(tree, b, radius);
    std::cout << num << std::endl;

    alglib::real_2d_array ans;
    alglib::real_1d_array distances;
    alglib::kdtreequeryresultsxy(tree, ans);
    alglib::kdtreequeryresultsdistances(tree, distances);

    ivec ret;
    ret.resize(num);
    for(int i = 0; i < num; i++) {
        std::cout << ans[i][3] << " " << distances[i] << std::endl;
        ret[i] = ans[i][3];
    }
    return ret;

}
*/
