#include "titanlib.h"

using namespace titanlib;

titanlib::KDTree::KDTree(vec lats, vec lons, CoordinateType type) {
    mLats = lats;
    mLons = lons;
    mType = type;
    vec x, y, z;

    titanlib::KDTree::convert_coordinates(mLats, mLons, x, y, z);
    for(int i = 0; i < mLats.size(); i++) {
        point p(x[i], y[i], z[i]);
        mTree.insert(std::make_pair(p, i));
    }
}

int titanlib::KDTree::get_num_neighbours(float lat, float lon, float radius, bool include_match) const {
    ivec indices = get_neighbours(lat, lon, radius, include_match);
    return indices.size();
}

ivec titanlib::KDTree::get_neighbours_with_distance(float lat, float lon, float radius, vec& distances, bool include_match) const {
    float x, y, z;
    titanlib::KDTree::convert_coordinates(lat, lon, x, y, z);
    ivec indices = get_neighbours(lat, lon, radius, include_match);

    int num = indices.size();
    distances.resize(num);
    for(int i = 0; i < num; i++) {
        float x1, y1, z1;
        titanlib::KDTree::convert_coordinates(mLats[indices[i]], mLons[indices[i]], x1, y1, z1);
        distances[i] = titanlib::KDTree::calc_distance(x, y, z, x1, y1, z1);
    }

    return indices;
}

ivec titanlib::KDTree::get_neighbours(float lat, float lon, float radius, bool include_match) const {
    float x, y, z;
    titanlib::KDTree::convert_coordinates(lat, lon, x, y, z);

    std::vector<value> results;
#if 1
    point p(x, y, z);
    within_radius r(p, radius);

    if(!include_match) {
        // Include match search
        is_not_equal s(p);
        mTree.query(boost::geometry::index::satisfies(r) && boost::geometry::index::satisfies(s), std::back_inserter(results));
    }
    else {
        mTree.query(boost::geometry::index::satisfies(r), std::back_inserter(results));
    }
    int num = results.size();

    ivec ret;
    ret.reserve(num);
    for(int i = 0; i < num; i++) {
        ret.push_back(results[i].second);
    }
#else
    box bx(point(x - radius, y - radius, z - radius), point(x + radius, y + radius, z + radius));
    mTree.query(boost::geometry::index::within(bx), std::back_inserter(results));
    int num = results.size();

    ivec ret;
    ret.reserve(num);
    for(int i = 0; i < num; i++) {
        float x1, y1, z1;
        titanlib::KDTree::convert_coordinates(mLats[results[i].second], mLons[results[i].second], x1, y1, z1);
        float dist = titanlib::KDTree::calc_distance(x, y, z, x1, y1, z1);
        if(dist <= radius) {
            ret.push_back(results[i].second);
        }
    }
#endif
    return ret;
}

ivec titanlib::KDTree::get_closest_neighbours(float lat, float lon, int num, bool include_match) const {
    float x, y, z;
    titanlib::KDTree::convert_coordinates(lat, lon, x, y, z);
    point p(x, y, z);

    std::vector<value> results;
    if(!include_match) {
        is_not_equal s(p);
        mTree.query(boost::geometry::index::nearest(p, num) && boost::geometry::index::satisfies(s), std::back_inserter(results));
    }
    else {
        mTree.query(boost::geometry::index::nearest(p, num), std::back_inserter(results));
    }
    int num_found = results.size();

    ivec ret;
    ret.reserve(num);
    for(int i = 0; i < num_found; i++) {
        ret.push_back(results[i].second);
    }
    return ret;
}
int titanlib::KDTree::get_nearest_neighbour(float lat, float lon, bool include_match) const {
    return get_closest_neighbours(lat, lon, 1, include_match)[0];
}
bool titanlib::KDTree::convert_coordinates(const vec& lats, const vec& lons, vec& x_coords, vec& y_coords, vec& z_coords) const {
    int N = lats.size();
    x_coords.resize(N);
    y_coords.resize(N);
    z_coords.resize(N);
    for(int i = 0; i < N; i++) {
        convert_coordinates(lats[i], lons[i], x_coords[i], y_coords[i], z_coords[i]);
    }

    return true;
}

bool titanlib::KDTree::convert_coordinates(float lat, float lon, float& x_coord, float& y_coord, float& z_coord) const {
    if(mType == titanlib::Cartesian) {
        x_coord = lon;
        y_coord = lat;
        z_coord = 0;
    }
    else {
        double lonr = M_PI / 180 * lon;
        double latr = M_PI / 180 * lat;
        // std::cout << lon << " " << lat << std::endl;
        x_coord = std::cos(latr) * std::cos(lonr) * titanlib::radius_earth;
        y_coord = std::cos(latr) * std::sin(lonr) * titanlib::radius_earth;
        z_coord = std::sin(latr) * titanlib::radius_earth;
    }
    return true;
}
float titanlib::KDTree::calc_distance(float lat1, float lon1, float lat2, float lon2, CoordinateType type) {
    if(type == titanlib::Cartesian) {
        float dx = lon1 - lon2;
        float dy = lat1 - lat2;
        return sqrt(dx * dx + dy * dy);
    }
    else if(type == titanlib::Geodetic){
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
}
float titanlib::KDTree::calc_distance_fast(float lat1, float lon1, float lat2, float lon2, CoordinateType type) {
    if(type == titanlib::Cartesian) {
        float dx = lon1 - lon2;
        float dy = lat1 - lat2;
        return sqrt(dx * dx + dy * dy);
    }
    else if(type == titanlib::Geodetic){
        double lat1r = deg2rad(lat1);
        double lat2r = deg2rad(lat2);
        double lon1r = deg2rad(lon1);
        double lon2r = deg2rad(lon2);
        float dx2 = pow(cos((lat1r+lat2r)/2),2)*(lon1r-lon2r)*(lon1r-lon2r);
        float dy2 = (lat1r-lat2r)*(lat1r-lat2r);
        return titanlib::radius_earth*sqrt(dx2+dy2);
    }
    else {
        throw std::runtime_error("Unknown coordinate type");
    }
}
float titanlib::KDTree::calc_distance_fast(const Point& p1, const Point& p2) {
    assert(p1.type == p2.type);
    return calc_distance_fast(p1.lat, p1.lon, p2.lat, p2.lon, p1.type);
}
float titanlib::KDTree::calc_distance(float x0, float y0, float z0, float x1, float y1, float z1) {
    return sqrt((x0 - x1)*(x0 - x1) + (y0 - y1)*(y0 - y1) + (z0 - z1)*(z0 - z1));
}
float titanlib::KDTree::deg2rad(float deg) {
   return (deg * M_PI / 180);
}
float titanlib::KDTree::rad2deg(float rad) {
   return (rad * 180 / M_PI);
}
vec titanlib::KDTree::get_lats() const {
    return mLats;
}
vec titanlib::KDTree::get_lons() const {
    return mLons;
}
int titanlib::KDTree::size() const {
    return mLats.size();
}
CoordinateType titanlib::KDTree::get_coordinate_type() const {
    return mType;
}
titanlib::KDTree& titanlib::KDTree::operator=(titanlib::KDTree other) {
    std::swap(mLats, other.mLats);
    std::swap(mLons, other.mLons);
    std::swap(mTree, other.mTree);
    std::swap(mType, other.mType);
    return *this;
}
titanlib::KDTree::KDTree(const titanlib::KDTree& other) {
    mLats = other.mLats;
    mLons = other.mLons;
    mTree = other.mTree;
    mType = other.mType;
}
titanlib::KDTree::within_radius::within_radius(point p, float radius)  {
    this->p = p;
    this->radius = radius;
}

bool titanlib::KDTree::within_radius::operator()(value const& v) const {
    float x0 = v.first.get<0>();
    float y0 = v.first.get<1>();
    float z0 = v.first.get<2>();
    float x1 = p.get<0>();
    float y1 = p.get<1>();
    float z1 = p.get<2>();
    return titanlib::KDTree::calc_distance(x0, y0, z0, x1, y1, z1) <= radius;
}
titanlib::KDTree::is_not_equal::is_not_equal(point p)  {
    this->p = p;
}

bool titanlib::KDTree::is_not_equal::operator()(value const& v) const {
    float x0 = v.first.get<0>();
    float y0 = v.first.get<1>();
    float z0 = v.first.get<2>();
    return p.get<0>() != x0 || p.get<1>() != y0 || p.get<2>() != z0;
}
