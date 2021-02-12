#include <cmath>
#include "titanlib.h"
#include <proj_api.h>
#include <math.h>

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

float titanlib::util::findKclosest(int k, const vec& array) {
    int n = array.size();
    if(n == 0) {
        throw std::runtime_error("Cannot compute quantile on empty array");
    }
    vec array_copy(n);
    // make a copy of the vector
    for(int i = 0; i < n; i++)
        array_copy[i] = array[i];
    std::sort(array_copy.begin(), array_copy.end());
    if(k > n) {
      k = n-1;
    } else {
      k = k-1;
    }
    float value = array_copy[k];
    return value;
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

float* titanlib::util::test_array(float* v, int n) {
    int count = 0;
    for(int i = 0; i < n; i++)
        count++;
    return v;
 }

//+ Set indices when inner and outer circles are present
bool titanlib::util::set_indices( const ivec& indices_global_outer_guess, 
                  const ivec& obs_test,
                  const ivec& dqcflags, 
                  const vec& dist_outer_guess, 
                  float inner_radius, 
                  int test_just_this,
                  ivec& indices_global_outer, 
                  ivec& indices_global_test, 
                  ivec& indices_outer_inner, 
                  ivec& indices_outer_test, 
                  ivec& indices_inner_test) {
/*-----------------------------------------------------------------------------
 Set all the vectors of indices linking the vectors operating over the
 different regions that we have defined: global, outer circle, inner circle,
 observations to test within the inner circle
 If the index test_just_this is set to a value greater than 0 (i.e. it 
 does point to an observation of the global vectors), then the only
 observation to test is test_just_this.
 
 Rules for an observation to be assigned to:
 1) outer circle. Within outer_radius from the centroid and flag!=1
 2) inner circle. In the outer circle AND within inner_radius from the centroid
 3) test. In the inner circle AND flag!=0 AND the user want to test it 

----------------------------------------------------------------------------*/
    int p_outer_guess = indices_global_outer_guess.size();
    indices_global_outer.reserve(p_outer_guess);
    indices_global_test.reserve(p_outer_guess);
    indices_outer_inner.reserve(p_outer_guess);
    indices_outer_test.reserve(p_outer_guess);
    indices_inner_test.reserve(p_outer_guess);
    int i = -1; // counter for indices on outer
    int l = -1; // counter for indices on inner
    for(int count=0; count<p_outer_guess; count++) {
        int g = indices_global_outer_guess[count];
        if(dqcflags[g] != 1 && g != test_just_this) {
            indices_global_outer.push_back(g);
            i++;
            if(dist_outer_guess[count] <= inner_radius) {
                indices_outer_inner.push_back(i);
                l++;
                if (test_just_this < 0 && dqcflags[g] != 0 && obs_test[g] == 1) {
                    indices_global_test.push_back(g);
                    indices_outer_test.push_back(i);
                    indices_inner_test.push_back(l);
                }
            }
        }
    }
    if ( test_just_this >= 0) {
        // add global index to global_outer indices as its last element
        indices_global_outer.push_back(test_just_this);
        // add global index to global_test as its only element
        indices_global_test.push_back(test_just_this);
        // adjust outer_inner, last element of outer_inner points to last element of outer 
        indices_outer_inner.push_back(indices_global_outer.size()-1);
        // adjust outer_test, only element of test points to last element of outer 
        indices_outer_test.push_back(indices_global_outer.size()-1);
        // adjust inner_test, only element of test points to last element of inner 
        indices_inner_test.push_back(indices_outer_inner.size()-1);
    }
    //
    return true;
}

