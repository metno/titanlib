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
    alglib::real_1d_array b;
    float x, y, z;
    titanlib::util::convert_coordinates(lat, lon, x, y, z);
    b.setlength(3);
    b[0] = x;
    b[1] = y;
    b[2] = z;
    int num = alglib::kdtreequeryrnn(mTree, b, radius);
    return num;
}
ivec titanlib::KDTree::get_neighbours(float lat, float lon, float radius) {
    alglib::real_1d_array b;
    float x, y, z;
    titanlib::util::convert_coordinates(lat, lon, x, y, z);
    b.setlength(3);
    b[0] = x;
    b[1] = y;
    b[2] = z;
    int num = alglib::kdtreequeryrnn(mTree, b, radius);

    ivec ret;
    alglib::real_2d_array ans;
    alglib::kdtreequeryresultsxy(mTree, ans);
    ret.resize(num);
    for(int i = 0; i < num; i++) {
        // std::cout << ans[i][3] << " " << distances[i] << std::endl;
        ret[i] = ans[i][3];
    }
    return ret;

}
