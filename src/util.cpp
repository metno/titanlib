#include <cmath>
#include "titanlib.h"
#include <proj_api.h>
#include <math.h>

using namespace titanlib;

bool titanlib::util::convert_coordinates(const vec& lats, const vec& lons, vec& x_coords, vec& y_coords, vec& z_coords) {
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
void titanlib::util::convert_to_proj(const vec& lats, const vec& lons, std::string proj4, vec& x_coords, vec& y_coords) {
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
vec titanlib::util::interpolate_to_points(const vec2& input_lats, const vec2& input_lons, const vec2& input_values, const vec& output_lats, const vec& output_lons) {
    assert(input_lats.size() > 0);
    assert(input_lats[0].size() > 0);
    vec output_values(output_lats.size(), 0);
    vec input_lats_flat(input_lats.size() * input_lats[0].size(), 0);
    vec input_lons_flat(input_lats.size() * input_lats[0].size(), 0);
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

float titanlib::util::compute_quantile(double quantile, const vec& array) {
    int n = array.size();
    if(n == 0) {
        throw std::runtime_error("Cannot compute quantile on empty array");
    }
    vec array_copy(n);
    // make a copy of the vector
    for(int i = 0; i < n; i++)
        array_copy[i] = array[i];
    float exact_q;
    std::sort(array_copy.begin(), array_copy.end());
    // get the quantile from sorted array
    int lowerIndex = floor(quantile * (n-1));
    int upperIndex = ceil(quantile * (n-1));
    float lowerValue = array_copy[lowerIndex];
    float upperValue = array_copy[upperIndex];
    float lowerQuantile = (float) lowerIndex / (n-1);
    float upperQuantile = (float) upperIndex / (n-1);
    if(lowerIndex == upperIndex) {
        exact_q = lowerValue;
    }
    else {
        assert(upperQuantile > lowerQuantile);
        assert(quantile >= lowerQuantile);
        float f = (quantile - lowerQuantile)/(upperQuantile - lowerQuantile);
        assert(f >= 0);
        assert(f <= 1);
        exact_q = lowerValue + (upperValue - lowerValue) * f;
    }

    return exact_q;
}

vec titanlib::util::subset(const vec& input, const ivec& indices) {
    vec output(indices.size());
    int size = indices.size();
    for(int i=0; i < size; i++) {
        int index = indices[i];
        assert(index < input.size());
        output[i] = input[index];
    }
    return output;
}

Points titanlib::util::subset(const Points& input, const ivec& indices) {
    int size = indices.size();
    vec ilats = input.get_lats();
    vec ilons = input.get_lons();
    vec ielevs = input.get_elevs();
    vec ilafs = input.get_lafs();
    vec lats(size);
    vec lons(size);
    vec elevs(size);
    vec lafs(size);
    for(int i=0; i < size; i++) {
        int index = indices[i];
        assert(index < input.size());
        lats[i] = ilats[index];
        lons[i] = ilons[index];
        elevs[i] = ielevs[index];
        lafs[i] = ilafs[index];
    }
    return Points(lats, lons, elevs, lafs);
}

float* titanlib::util::test_array(float* v, int n) {
    int count = 0;
    for(int i = 0; i < n; i++)
        count++;
    return v;
 }
