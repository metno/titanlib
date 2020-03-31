#include <cmath>
#include "titanlib.h"
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
fvec titanlib::util::interpolate_to_points(const fvec2& input_lats, const fvec2& input_lons, const fvec2& input_values, const fvec& output_lats, const fvec& output_lons) {
    assert(input_lats.size() > 0);
    assert(input_lats[0].size() > 0);
    fvec output_values(output_lats.size(), 0);
    fvec input_lats_flat(input_lats.size() * input_lats[0].size(), 0);
    fvec input_lons_flat(input_lats.size() * input_lats[0].size(), 0);
    int X = input_lats[0].size();
    int Y = input_lats.size();
    int count = 0;
    for(int i = 0; i < Y; i++) {
        for(int j = 0; j < X; j++) {
            input_lats_flat[count] = input_lats[i][j];
            input_lons_flat[count] = input_lons[i][j];
            count++;
        }
    }
    titanlib::KDTree tree(input_lats_flat, input_lons_flat);
    for(int i = 0; i < output_lats.size(); i++) {
        int index = tree.get_nearest_neighbour(output_lats[i], output_lons[i], true);
        int index_x = index % X;
        int index_y = index / X;
        output_values[i] = input_values[index_y][index_x];
    }

    return output_values;
}
