#include "titanlib.h"

titanlib::KDTree::KDTree(const fvec& lats, const fvec& lons) {
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
    }
    alglib::kdtreebuild(a, nx, ny, normtype, mTree);
}

int titanlib::KDTree::get_num_neighbours(float lat, float lon, float radius) {
    alglib::real_1d_array b = titanlib::KDTree::ll2ar(lat, lon);
    int num = alglib::kdtreequeryrnn(mTree, b, radius, false);
    return num;
}

ivec titanlib::KDTree::get_neighbours(float lat, float lon, float radius) {
    alglib::real_1d_array b = titanlib::KDTree::ll2ar(lat, lon);
    int num = alglib::kdtreequeryrnn(mTree, b, radius, false);

    ivec ret;
    alglib::real_2d_array ans;
    alglib::kdtreequeryresultsxy(mTree, ans);
    ret.resize(num);
    for(int i = 0; i < num; i++) {
        ret[i] = ans[i][3];
    }
    return ret;

}

ivec titanlib::KDTree::get_neighbours_with_distance(float lat, float lon, float radius, fvec& distances) {
    alglib::real_1d_array b = titanlib::KDTree::ll2ar(lat, lon);
    int num = alglib::kdtreequeryrnn(mTree, b, radius, false);

    ivec ret;
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
    return ret;
}

ivec titanlib::KDTree::get_closest_neighbours(float lat, float lon, int num) {

    alglib::real_1d_array b = titanlib::KDTree::ll2ar(lat, lon);
    int num_found = alglib::kdtreequeryknn(mTree, b, num, false);

    ivec ret;
    alglib::real_2d_array ans;
    alglib::kdtreequeryresultsxy(mTree, ans);
    ret.resize(num_found);
    for(int i = 0; i < num_found; i++) {
        ret[i] = ans[i][3];
    }

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

    ivec ret;
    alglib::real_2d_array ans;
    alglib::kdtreequeryresultsxy(mTree, ans);
    return ans[0][3];
}
