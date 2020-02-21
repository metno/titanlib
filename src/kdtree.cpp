#include "titanlib.h"
#define USE_ALGLIB 0

titanlib::KDTree::KDTree(const fvec& lats, const fvec& lons) {
    mLats = lats;
    mLons = lons;
    fvec x, y, z;
    int nx = 3;
    int ny = 1;
    int normtype = 2;
    titanlib::util::convert_coordinates(lats, lons, x, y, z);

    alglib::real_2d_array a;
    a.setlength(lats.size(), 4);
    for(int i = 0; i < lats.size(); i++) {
        a[i][0] = x[i];
        a[i][1] = y[i];
        a[i][2] = z[i];
        a[i][3] = i;
        point p(x[i], y[i], z[i]);
        mTree2.insert(std::make_pair(p, i));

    }
    alglib::kdtreebuild(a, nx, ny, normtype, mTree);
}

int titanlib::KDTree::get_num_neighbours(float lat, float lon, float radius) {
    alglib::real_1d_array b = titanlib::KDTree::ll2ar(lat, lon);
    int num = alglib::kdtreequeryrnn(mTree, b, radius, false);
    return num;
}

ivec titanlib::KDTree::get_neighbours(float lat, float lon, float radius) {
    ivec ret;
#if USE_ALGLIB
    alglib::real_1d_array b = titanlib::KDTree::ll2ar(lat, lon);
    int num = alglib::kdtreequeryrnn(mTree, b, radius, false);

    alglib::real_2d_array ans;
    alglib::kdtreequeryresultsxy(mTree, ans);
    ret.resize(num);
    for(int i = 0; i < num; i++) {
        ret[i] = ans[i][3];
    }
#else
    std::vector<value> result_n;
    float x, y, z;
    titanlib::util::convert_coordinates(lat, lon, x, y, z);
    box bx(point(x - radius, y - radius, z - radius), point(x + radius, y + radius, z + radius));
    mTree2.query(boost::geometry::index::within(bx), std::back_inserter(result_n));
    int num = result_n.size();

    ivec ret2;
    ret2.resize(num);
    for(int i = 0; i < num; i++) {
        ret2[i] = result_n[i].second;
    }
#endif
    return ret2;

}

ivec titanlib::KDTree::get_neighbours_with_distance(float lat, float lon, float radius, fvec& distances) {
    ivec ret;
#if USE_ALGLIB
    int
    alglib::real_1d_array b = titanlib::KDTree::ll2ar(lat, lon);
    int num = alglib::kdtreequeryrnn(mTree, b, radius, false);

    alglib::real_2d_array ans;
    alglib::kdtreequeryresultsxy(mTree, ans);
    ret.resize(num);
    for(int i = 0; i < num; i++) {
        ret[i] = ans[i][3];
    }

    alglib::real_1d_array rdist;
    alglib::kdtreequeryresultsdistances(mTree, rdist);
    distances.resize(num);
    for(int i = 0; i < num; i++) {
        distances[i] = rdist[i];
    }
#else
    std::vector<value> result_n;
    float x, y, z;
    titanlib::util::convert_coordinates(lat, lon, x, y, z);
    box bx(point(x - radius, y - radius, z - radius), point(x + radius, y + radius, z + radius));
    mTree2.query(boost::geometry::index::within(bx), std::back_inserter(result_n));
    int num = result_n.size();

    distances.resize(num);
    for(int i = 0; i < num; i++) {
        distances[i] = titanlib::util::calc_distance(lat, lon, mLats[result_n[i].second], mLons[result_n[i].second]);
    }

    ivec ret2;
    ret2.resize(num);
    for(int i = 0; i < num; i++) {
        ret2[i] = result_n[i].second;
    }

#endif
    return ret;
}

ivec titanlib::KDTree::get_closest_neighbours(float lat, float lon, int num) {
    ivec ret;
#if USE_ALGLIB
    alglib::real_1d_array b = titanlib::KDTree::ll2ar(lat, lon);
    int num_found = alglib::kdtreequeryknn(mTree, b, num, false);

    alglib::real_2d_array ans;
    alglib::kdtreequeryresultsxy(mTree, ans);
    ret.resize(num_found);
    for(int i = 0; i < num_found; i++) {
        ret[i] = ans[i][3];
    }
#else
    std::vector<value> result_n;
    float x, y, z;
    titanlib::util::convert_coordinates(lat, lon, x, y, z);
    mTree2.query(boost::geometry::index::nearest(point(x, y, z), num), std::back_inserter(result_n));
    int num_found = result_n.size();

    ret.resize(num_found);
    for(int i = 0; i < num_found; i++) {
        ret[i] = result_n[i].second;
    }
#endif
    return ret;
}
alglib::real_1d_array titanlib::KDTree::ll2ar(float lat, float lon) {
    alglib::real_1d_array b;
    float x, y, z;
    titanlib::util::convert_coordinates(lat, lon, x, y, z);
    b.setlength(3);
    b[0] = x;
    b[1] = y;
    b[2] = z;
    return b;
}
int titanlib::KDTree::get_nearest_neighbour(float lat, float lon) {
    alglib::real_1d_array b = titanlib::KDTree::ll2ar(lat, lon);
    int num_found = alglib::kdtreequeryknn(mTree, b, 1, false);
    assert(num_found == 1);

    ivec ret;
    alglib::real_2d_array ans;
    alglib::kdtreequeryresultsxy(mTree, ans);
    return ans[0][3];
}
