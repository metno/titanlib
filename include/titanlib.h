#ifndef TITANLIB_H
#define TITANLIB_H
#include <iostream>
#include <vector>
#include <assert.h>
#include <boost/geometry.hpp>
#include <boost/geometry/geometries/point.hpp>
#include <boost/geometry/geometries/box.hpp>
#include <boost/geometry/index/rtree.hpp>

#define TITANLIB_VERSION "0.1.0"
#define __version__ TITANLIB_VERSION

typedef std::vector<int> ivec;
typedef std::vector<float> fvec;
typedef std::vector<double> dvec;
typedef std::vector<fvec> fvec2;

/** Titanlib
*/
namespace titanlib {
    /**
     * @return Titanlib version
     */
    std::string version();

    /** Spatial Consistency Test
        @param lats vector of latitudes
        @param flags output vector of flags
     */
    ivec sct(const fvec& lats,
            const fvec& lons,
            const fvec& elevs,
            const fvec& values,
            int nmin,
            int nmax,
            int nminprof,
            float dzmin,
            float dhmin,
            float dz,
            const fvec& t2pos,
            const fvec& t2neg,
            const fvec& eps2,
            fvec& sct,
            fvec& rep);

    /** Spatial Consistency Test
     *  @param nminprof Minimum number of observations to compute vertical profile
     *  @param radius Select observations within this radius [m]
     *  @param dzmin Minimum elevation difference to compute vertical profile [m]
     *  @param dhmin Minimum horizontal decorrelation length [m]
     *  @param dz Vertical decorrelation length [m]
     *  @param pos Positive deviation allowed
     *  @param neg Negative deviation allowed
     *  @param eps2
     *  @param sct
     *  @param flags
     */
    ivec sct_rewrite(const fvec& lats,
            const fvec& lons,
            const fvec& elevs,
            const fvec& values,
            int minnumobs,
            int maxnumobs,
            double inner_radius,
            double outer_radius,
            int nminprof,
            double dzmin,
            double dhmin,
            float dz,
            const fvec& t2pos,
            const fvec& t2neg,
            const fvec& eps2);

    /** Range check. Checks observation is within the ranges given
     *  @param values vector of observations
     *  @param min min allowed value
     *  @param max max allowed value
     *  @param flags vector of return flags
     */
    ivec range_check(const fvec& values,
            const fvec& min,
            const fvec& max);

    ivec range_check_climatology(const fvec& lats,
            const fvec& lons,
            const fvec& elevs,
            const fvec& values,
            int unixtime,
            const fvec& plus,
            const fvec& minus);

    /** Buddy check. Compares a station to all its neighbours within a certain distance
     *  @param lats vector of latitudes [deg]
     *  @param lons vector of longitudes [deg]
     *  @param elevs vector of elevations [m]
     *  @param values vector of observation values
     *  @param radius search radius [m]
     *  @param buddies_min the minimum number of buddies a station can have
     *  @param thresholds the threshold for flagging a station
     *  @param diff_elev_max the maximum difference in elevation for a buddy (if negative will not check for heigh difference)
     *  @param elev_gradient the multiplier used to adjust the values based on elevation difference from the comparison value
     *  @param min_std 
     *  @param num_iterations 
     *  @param obs_to_check the observations that will be checked (since can pass in observations that will not be checked)
     *  @param flags vector of return flags
     */
    ivec buddy_check(const fvec& lats,
            const fvec& lons,
            const fvec& elevs,
            const fvec& values,
            const fvec& radius,
            const ivec& buddies_min,
            const fvec& thresholds,
            float diff_elev_max,
            float elev_gradient,
            float min_std,
            int num_iterations,
            const ivec obs_to_check = ivec());

    ivec buddy_event_check(const fvec& lats,
            const fvec& lons,
            const fvec& elevs,
            const fvec& values,
            const fvec& radius,
            const ivec& buddies_min,
            const fvec& event_thresholds,
            const fvec& thresholds,
            float diff_elev_max,
            float elev_gradient,
            const ivec obs_to_check = ivec());

    ivec first_guess_check(
        const fvec& values,
        const fvec& fg_values,
        const fvec& fg_neg,
        const fvec& fg_pos);

    /** Isolation check. Checks that a station is not located alone
     *  @param lats vector of latitudes [deg]
     *  @param lons vector of longitudes [deg]
     *  @param nmin required number of observations
     *  @param radius search radius [m]
     *  @param flags vector of return flags
     */
    ivec isolation_check(const fvec& lats,
            const fvec& lons,
            int nmin,
            float radius);

    /** Isolation check with elevation. Checks that a station is not located alone
     *  @param lats vector of latitudes [deg]
     *  @param lons vector of longitudes [deg]
     *  @param elevs vector of elevations [m]
     *  @param nmin required number of observations
     *  @param radius search radius [m]
     *  @param dz vertical search radius [m]
     *  @param flags vector of return flags
     */
    ivec isolation_check(const fvec& lats,
            const fvec& lons,
            const fvec& elevs,
            int nmin,
            float radius,
            float dz);

    namespace util {
        /** Convert lat/lons to 3D cartesian coordinates with the centre of the earth as the origin
         *  @param lats vector of latitudes [deg]
         *  @param lons vector of longitudes [deg]
         *  @param x_coords vector of x-coordinates [m]
         *  @param y_coords vector of y-coordinates [m]
         *  @param z_coords vector of z-coordinates [m]
         */
        bool convert_coordinates(const fvec& lats, const fvec& lons, fvec& x_coords, fvec& y_coords, fvec& z_coords);

        /** Same as above, but convert a single lat/lon to 3D cartesian coordinates
         *  @param lat latitude [deg]
         *  @param lon longitude [deg]
         *  @param x_coord x-coordinate [m]
         *  @param y_coord y-coordinate [m]
         *  @param z_coord z-coordinate [m]
         */
        bool convert_coordinates(float lat, float lon, float& x_coord, float& y_coord, float& z_coord);

        void convert_to_proj(const fvec& lats, const fvec& lons, std::string proj4, fvec& x_coords, fvec& y_coords);
        fvec interpolate_to_points(const fvec2& input_lats, const fvec2& input_lons, const fvec2& input_values, const fvec& output_lats, const fvec& output_lons);
        float deg2rad(float deg);
        float calc_distance(float lat1, float lon1, float lat2, float lon2);
        float calc_distance(float x0, float y0, float z0, float x1, float y1, float z1);
    }

    // ivec nearest_neighbours(const fvec& lats, const fvec& lons, float radius, float lat, float lon);

    // bool prioritize(const fvec& values, const ivec& priority, float distance, ivec& flags);

    /** Represents point and their observed values */
    class Dataset {
        public:
            Dataset(fvec ilats, fvec ilons, fvec ielevs, fvec ivalues);
            /** Perform the range check on the dataset
             *  @param indices Only perform the test on these indices
             */
            bool range_check(const fvec& min, const fvec& max, const ivec indices=ivec());
            bool range_check_climatology(int unixtime, const fvec& plus, const fvec& minus, const ivec indices=ivec());
            bool sct(int nmin, int nmax, int nminprof, float dzmin, float dhmin, float dz, const fvec& t2pos, const fvec& t2neg, const fvec& eps2, fvec& sct, fvec& rep, const ivec indices=ivec());
            bool buddy_check(const fvec& radius, const ivec& buddies_min, const fvec& thresholds, float diff_elev_max, float elev_gradient, float min_std, int num_iteratiowns, const ivec& obs_to_check, const ivec indices=ivec());

            fvec lats;
            fvec lons;
            fvec elevs;
            fvec values;
            ivec flags;
        private:
            template <class T> T subset(const T& array, const ivec& indices) {
                if(array.size() == 1) {
                    T new_array = array;
                    return new_array;
                }
                else {
                    T new_array(indices.size());
                    for(int i = 0; i < indices.size(); i++) {
                        new_array[i] = array[indices[i]];
                    }
                    return new_array;
                }
            }
            template <class T> void unsubset(const T& array, T& orig_array, const ivec& indices) {
                assert(array.size() == indices.size());
                for(int i = 0; i < indices.size(); i++) {
                    orig_array[indices[i]] = array[i];
                }
            }
    };

    class KDTree {
        public:
            KDTree(const fvec& lats, const fvec& lons);

            /** Find single nearest points
             *  @param lat Latitude of lookup-point
             *  @param lon Longitude of lookup-point
             */
            int get_nearest_neighbour(float lat, float lon, bool include_match);
            ivec get_nearest_neighbour(const fvec& lats, const fvec& lons, bool include_match);

            /** Find all points with a radius
             *  @param lat Latitude of lookup-point
             *  @param lon Longitude of lookup-point
             *  @param radius Lookup radius (use 0 for unlimited) [m]
             *  @param max_num Maximum number of points (use 0 for unlimited)
             *  @param include_match Should an exact match be included?
             */
            ivec get_neighbours(float lat, float lon, float radius, int max_num, bool include_match);

            /** Find all points with a radius
             *  @param lat Latitude of lookup-point
             *  @param lon Longitude of lookup-point
             *  @param radius Lookup radius [m]
             *  @param distances Vector to store separation distances [m]
             */
            ivec get_neighbours_with_distance(float lat, float lon, float radius, int max_num, bool include_match, fvec& distances);

            /** Find the number of points within a radius
             *  @param lat Latitude of lookup-point
             *  @param lon Longitude of lookup-point
             *  @param radius Lookup radius [m]
             */
            int get_num_neighbours(float lat, float lon, float radius, int max_num, bool include_match);
        private:
            typedef boost::geometry::model::point<float, 3, boost::geometry::cs::cartesian> point;
            typedef std::pair<point, unsigned> value;
            typedef boost::geometry::model::box<point> box;
            fvec mLats;
            fvec mLons;
            boost::geometry::index::rtree< value, boost::geometry::index::quadratic<16> > mTree;

    };
}
#endif
