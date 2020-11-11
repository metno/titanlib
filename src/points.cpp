#include "titanlib.h"

using namespace titanlib;

titanlib::Points::Points() {
    // TODO: Deal with the empty case. Can't really find nearest neighbours then

}
titanlib::Points::Points(vec lats, vec lons, vec elevs, vec lafs, CoordinateType type) {
    int N = lats.size();
    if(lons.size() != N)
        throw std::invalid_argument("Cannot create points with unequal lat and lon sizes");
    if(elevs.size() != 0 && elevs.size() != N)
        throw std::invalid_argument("'elevs' must either be size 0 or the same size at lats/lons");
    if(lafs.size() != 0 && lafs.size() != N)
        throw std::invalid_argument("'lafs' must either be size 0 or the same size at lats/lons");
    mLats = lats;
    mLons = lons;
    mElevs = elevs;
    mLafs = lafs;
    KDTree tree = KDTree(lats, lons, type);
    mTree = tree;
    if(mElevs.size() != N) {
        mElevs.clear();
        mElevs.resize(N, titanlib::MV);
    }
    if(mLafs.size() != N) {
        mLafs.clear();
        mLafs.resize(N, titanlib::MV);
    }
}
titanlib::Points::Points(KDTree tree, vec elevs, vec lafs) {
    mElevs = elevs;
    mLafs = lafs;
    mTree = tree;
    mLats = tree.get_lats();
    mLons = tree.get_lons();
}

int titanlib::Points::get_num_neighbours(float lat, float lon, float radius, bool include_match) const {
    return mTree.get_num_neighbours(lat, lon, radius, include_match);
}

ivec titanlib::Points::get_neighbours_with_distance(float lat, float lon, float radius, vec& distances, bool include_match) const {
    return mTree.get_neighbours_with_distance(lat, lon, radius, distances, include_match);
}

ivec titanlib::Points::get_neighbours(float lat, float lon, float radius, bool include_match) const {
    return mTree.get_neighbours(lat, lon, radius, include_match);
}

ivec titanlib::Points::get_closest_neighbours(float lat, float lon, int num, bool include_match) const {
    return mTree.get_closest_neighbours(lat, lon, num, include_match);
}
int titanlib::Points::get_nearest_neighbour(float lat, float lon, bool include_match) const {
    ivec I = get_closest_neighbours(lat, lon, 1, include_match);
    if(I.size() > 0)
        return I[0];
    else
        return -1;
}
vec titanlib::Points::get_lats() const {
    return mLats;
}
vec titanlib::Points::get_lons() const {
    return mLons;
}
vec titanlib::Points::get_elevs() const {
    return mElevs;
}
vec titanlib::Points::get_lafs() const {
    return mLafs;
}
int titanlib::Points::size() const {
    return mLats.size();
}
titanlib::Points& titanlib::Points::operator=(titanlib::Points other) {
    std::swap(mLats, other.mLats);
    std::swap(mLons, other.mLons);
    std::swap(mElevs, other.mElevs);
    std::swap(mLafs, other.mLafs);
    std::swap(mTree, other.mTree);
    return *this;
}
titanlib::Points::Points(const titanlib::Points& other) {
    mLats = other.mLats;
    mLons = other.mLons;
    mElevs = other.mElevs;
    mLafs = other.mLafs;
    mTree = KDTree(mLats, mLons, mTree.get_coordinate_type());
}
CoordinateType titanlib::Points::get_coordinate_type() const {
    return mTree.get_coordinate_type();
}
