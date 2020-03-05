#include "titanlib.h"

titanlib::KDTree::KDTree(const fvec& lats, const fvec& lons) {
    mLats = lats;
    mLons = lons;
    fvec x, y, z;
    titanlib::util::convert_coordinates(lats, lons, x, y, z);

    for(int i = 0; i < lats.size(); i++) {
        point p(x[i], y[i], z[i]);
        mTree.insert(std::make_pair(p, i));
    }
}

int titanlib::KDTree::get_num_neighbours(float lat, float lon, float radius, int max_num, bool include_match) {
    ivec indices = get_neighbours(lat, lon, radius, max_num, include_match);
    return indices.size();
}

ivec titanlib::KDTree::get_neighbours_with_distance(float lat, float lon, float radius, int max_num, bool include_match, fvec& distances) {
    float x, y, z;
    titanlib::util::convert_coordinates(lat, lon, x, y, z);
    ivec indices = get_neighbours(lat, lon, radius, max_num, include_match);

    int num = indices.size();
    distances.resize(num);
    for(int i = 0; i < num; i++) {
        float x1, y1, z1;
        titanlib::util::convert_coordinates(mLats[indices[i]], mLons[indices[i]], x1, y1, z1);
        distances[i] = titanlib::util::calc_distance(x, y, z, x1, y1, z1);
    }

    return indices;
}

ivec titanlib::KDTree::get_neighbours(float lat, float lon, float radius, int max_num, bool include_match) {
    assert(max_num > 0 || radius > 0);

    // Nearest num search
    float x, y, z;
    titanlib::util::convert_coordinates(lat, lon, x, y, z);
    point p(x, y, z);

    // Radius search
    box bx(point(x - radius, y - radius, z - radius), point(x + radius, y + radius, z + radius));

    std::vector<value> results;
    if(!include_match) {
        // Include match search
        struct is_not_equal {
            bool operator()(value const& v) const {
                float x0 = v.first.get<0>();
                float y0 = v.first.get<1>();
                float z0 = v.first.get<2>();
                return p.get<0>() != x0 || p.get<1>() != y0 || p.get<2>() != z0;
            };
            point p;
        };
        is_not_equal s;
        s.p = p;

        if(max_num == 0 && radius > 0)
            mTree.query(boost::geometry::index::within(bx) && boost::geometry::index::satisfies(s), std::back_inserter(results));
        else if(max_num > 0 && radius == 0)
            mTree.query(boost::geometry::index::nearest(p, max_num) && boost::geometry::index::satisfies(s), std::back_inserter(results));
        else if(max_num > 0 && radius > 0)
            mTree.query(boost::geometry::index::nearest(p, max_num) && boost::geometry::index::within(bx) && boost::geometry::index::satisfies(s), std::back_inserter(results));
    }
    else {
        if(max_num == 0 && radius > 0)
            mTree.query(boost::geometry::index::within(bx), std::back_inserter(results));
        else if(max_num > 0 && radius == 0)
            mTree.query(boost::geometry::index::nearest(p, max_num), std::back_inserter(results));
        else if(max_num > 0 && radius > 0)
            mTree.query(boost::geometry::index::nearest(p, max_num) && boost::geometry::index::within(bx), std::back_inserter(results));
    }
    int num_found = results.size();

    ivec ret;
    ret.reserve(num_found);
    for(int i = 0; i < num_found; i++) {
        ret.push_back(results[i].second);
    }
    return ret;
}

int titanlib::KDTree::get_nearest_neighbour(float lat, float lon, bool include_match) {
    ivec neighbours = get_neighbours(lat, lon, 0, 1, include_match);
    assert(neighbours.size() == 1);
    return neighbours[0];
}
