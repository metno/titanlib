#include <cmath>
#include "titanlib.h"
#include <libalglib/interpolation.h>
#include <proj_api.h>
#include <math.h>

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
void titanlib::util::convert_to_proj(const fvec& lats, const fvec& lons, std::string proj4, fvec& x_coords, fvec& y_coords) {
    int N = lats.size();
    x_coords.resize(N);
    y_coords.resize(N);

    projPJ Pproj = pj_init_plus(proj4.c_str());
    projPJ Plonglat = pj_init_plus("+proj=longlat +ellps=clrk66");

    for(int i = 0; i < N; i++) {
        double lat0 = deg2rad(lats[i]);
        double lon0 = deg2rad(lons[i]);
        pj_transform(Plonglat, Pproj, 1, 1, &lon0, &lat0, NULL);
        x_coords[i] = lon0;
        y_coords[i] = lat0;
    }
    pj_free(Pproj);
    pj_free(Plonglat);
}
float titanlib::util::calc_distance(float lat1, float lon1, float lat2, float lon2) {
    if(!(fabs(lat1) <= 90 && fabs(lat2) <= 90 && fabs(lon1) <= 360 && fabs(lon2) <= 360)) {
        // std::cout << " Cannot calculate distance, invalid lat/lon: (" << lat1 << "," << lon1 << ") (" << lat2 << "," << lon2 << ")";
        // std::cout << '\n';
    }
    if(lat1 == lat2 && lon1 == lon2)
        return 0;

    double lat1r = deg2rad(lat1);
    double lat2r = deg2rad(lat2);
    double lon1r = deg2rad(lon1);
    double lon2r = deg2rad(lon2);
    double radiusEarth = 6.378137e6;

    double ratio = cos(lat1r)*cos(lon1r)*cos(lat2r)*cos(lon2r)
                   + cos(lat1r)*sin(lon1r)*cos(lat2r)*sin(lon2r)
                   + sin(lat1r)*sin(lat2r);
    double dist = acos(ratio)*radiusEarth;
    return (float) dist;
}
float titanlib::util::calc_distance(float x0, float y0, float z0, float x1, float y1, float z1) {
    return sqrt((x0 - x1)*(x0 - x1) + (y0 - y1)*(y0 - y1) + (z0 - z1)*(z0 - z1));
}
float titanlib::util::deg2rad(float deg) {
   return (deg * M_PI / 180);
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
