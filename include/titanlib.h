#ifndef TITANLIB_H
#define TITANLIB_H
#include <iostream>
#include <vector>
#include <assert.h>
#include <boost/geometry.hpp>
#include <boost/geometry/geometries/point.hpp>
#include <boost/geometry/geometries/box.hpp>
#include <boost/geometry/index/rtree.hpp>
#ifdef _OPENMP
    #include <omp.h>
#endif

#define TITANLIB_VERSION "0.1.0"
#define __version__ TITANLIB_VERSION

typedef std::vector<int> ivec;
typedef std::vector<float> vec;
typedef std::vector<double> dvec;
typedef std::vector<vec> vec2;

/** Titanlib
*/
namespace titanlib {
    /**
     * @return Titanlib version
     */
    std::string version();

    float* test_array(float* v, int n);

    /** Spatial Consistency Test
        @param lats vector of latitudes
        @param flags output vector of flags
     */
    ivec sct_old(const vec& lats,
            const vec& lons,
            const vec& elevs,
            const vec& values,
            int nmin,
            int nmax,
            int num_min_prof,
            float dzmin,
            float dhmin,
            float dz,
            const vec& pos,
            const vec& neg,
            const vec& eps2,
            vec& sct,
            vec& rep,
            ivec& boxids);

    /** Spatial Consistency Test
     *  @param num_min_prof Minimum number of observations to compute vertical profile
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
    ivec sct(const vec& lats,
            const vec& lons,
            const vec& elevs,
            const vec& values,
            int num_min,
            int num_max,
            double inner_radius,
            double outer_radius,
            int num_iterations,
            int num_min_prof,
            double dzmin,
            double dhmin,
            float dz,
            const vec& pos,
            const vec& neg,
            const vec& eps2,
            vec& sct,
            vec& rep);

    /** Range check. Checks observation is within the ranges given
     *  @param values vector of observations
     *  @param min min allowed value
     *  @param max max allowed value
     *  @param flags vector of return flags
     */
    ivec range_check(const vec& values,
            const vec& min,
            const vec& max);

    ivec range_check_climatology(const vec& lats,
            const vec& lons,
            const vec& elevs,
            const vec& values,
            int unixtime,
            const vec& pos,
            const vec& neg);

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
    ivec buddy_check(const vec& lats,
            const vec& lons,
            const vec& elevs,
            const vec& values,
            const vec& radius,
            const ivec& buddies_min,
            const vec& thresholds,
            float diff_elev_max,
            float elev_gradient,
            float min_std,
            int num_iterations,
            const ivec obs_to_check = ivec());

    ivec buddy_event_check(const vec& lats,
            const vec& lons,
            const vec& elevs,
            const vec& values,
            const vec& radius,
            const ivec& buddies_min,
            const vec& event_thresholds,
            const vec& thresholds,
            float diff_elev_max,
            float elev_gradient,
            const ivec obs_to_check = ivec());

    ivec first_guess_check(
        const vec& values,
        const vec& first_guess,
        const vec& neg,
        const vec& pos);

    /** Isolation check. Checks that a station is not located alone
     *  @param lats vector of latitudes [deg]
     *  @param lons vector of longitudes [deg]
     *  @param num_min required number of observations
     *  @param radius search radius [m]
     *  @param flags vector of return flags
     */
    ivec isolation_check(const vec& lats,
            const vec& lons,
            int num_min,
            float radius);

    /** Isolation check with elevation. Checks that a station is not located alone
     *  @param lats vector of latitudes [deg]
     *  @param lons vector of longitudes [deg]
     *  @param elevs vector of elevations [m]
     *  @param num_min required number of observations
     *  @param radius search radius [m]
     *  @param dz vertical search radius [m]
     *  @param flags vector of return flags
     */
    ivec isolation_check(const vec& lats,
            const vec& lons,
            const vec& elevs,
            int num_min,
            float radius,
            float dz);

    /** Set the number of OpenMP threads to use. Overwrides OMP_NUM_THREAD env variable. */
    void set_omp_threads(int num);

    /** Sets the number of OpenMP threads to 1 if OMP_NUM_THREADS undefined */
    void initialize_omp();

    namespace util {
        double clock();
        /** Convert lat/lons to 3D cartesian coordinates with the centre of the earth as the origin
         *  @param lats vector of latitudes [deg]
         *  @param lons vector of longitudes [deg]
         *  @param x_coords vector of x-coordinates [m]
         *  @param y_coords vector of y-coordinates [m]
         *  @param z_coords vector of z-coordinates [m]
         */
        bool convert_coordinates(const vec& lats, const vec& lons, vec& x_coords, vec& y_coords, vec& z_coords);

        /** Same as above, but convert a single lat/lon to 3D cartesian coordinates
         *  @param lat latitude [deg]
         *  @param lon longitude [deg]
         *  @param x_coord x-coordinate [m]
         *  @param y_coord y-coordinate [m]
         *  @param z_coord z-coordinate [m]
         */
        bool convert_coordinates(float lat, float lon, float& x_coord, float& y_coord, float& z_coord);

        void convert_to_proj(const vec& lats, const vec& lons, std::string proj4, vec& x_coords, vec& y_coords);
        vec interpolate_to_points(const vec2& input_lats, const vec2& input_lons, const vec2& input_values, const vec& output_lats, const vec& output_lons);
        float deg2rad(float deg);
        float calc_distance(float lat1, float lon1, float lat2, float lon2);
        float calc_distance(float x0, float y0, float z0, float x1, float y1, float z1);

        float compute_quantile(double quantile, const vec& array);
        vec subset(const vec& input, const ivec& indices);
    }

    // ivec nearest_neighbours(const vec& lats, const vec& lons, float radius, float lat, float lon);

    // bool prioritize(const vec& values, const ivec& priority, float distance, ivec& flags);

    /** Represents point and their observed values */
    class Dataset {
        public:
            Dataset(vec ilats, vec ilons, vec ielevs, vec ivalues);
            /** Perform the range check on the dataset
             *  @param indices Only perform the test on these indices
             */
            bool range_check(const vec& min, const vec& max, const ivec indices=ivec());
            bool range_check_climatology(int unixtime, const vec& pos, const vec& neg, const ivec indices=ivec());
            bool sct(int num_min, int num_max, double inner_radius, double outer_radius, int num_iterations, int num_min_prof, double dzmin, double dhmin, float dz, const vec& t2pos, const vec& t2neg, const vec& eps2, vec& sct, vec& rep, const ivec indices=ivec());
            bool buddy_check(const vec& radius, const ivec& buddies_min, const vec& thresholds, float diff_elev_max, float elev_gradient, float min_std, int num_iteratiowns, const ivec& obs_to_check, const ivec indices=ivec());
            bool isolation_check(int num_min, float radius, float dz);

            vec lats;
            vec lons;
            vec elevs;
            vec values;
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
            KDTree(const vec& lats, const vec& lons);

            /** Find single nearest points
             *  @param lat Latitude of lookup-point
             *  @param lon Longitude of lookup-point
             */
            int get_nearest_neighbour(float lat, float lon, bool include_match);
            ivec get_nearest_neighbour(const vec& lats, const vec& lons, bool include_match);

            /** Find all points with a radius
             *  @param lat Latitude of lookup-point
             *  @param lon Longitude of lookup-point
             *  @param radius Lookup radius (use 0 for unlimited) [m]
             *  @param num_max Maximum number of points (use 0 for unlimited)
             *  @param include_match Should an exact match be included?
             */
            ivec get_neighbours(float lat, float lon, float radius, int num_max, bool include_match);

            /** Find all points with a radius
             *  @param lat Latitude of lookup-point
             *  @param lon Longitude of lookup-point
             *  @param radius Lookup radius [m]
             *  @param distances Vector to store separation distances [m]
             */
            ivec get_neighbours_with_distance(float lat, float lon, float radius, int num_max, bool include_match, vec& distances);

            /** Find the number of points within a radius
             *  @param lat Latitude of lookup-point
             *  @param lon Longitude of lookup-point
             *  @param radius Lookup radius [m]
             */
            int get_num_neighbours(float lat, float lon, float radius, int num_max, bool include_match);
        private:
            typedef boost::geometry::model::point<float, 3, boost::geometry::cs::cartesian> point;
            typedef std::pair<point, unsigned> value;
            typedef boost::geometry::model::box<point> box;
            vec mLats;
            vec mLons;
            boost::geometry::index::rtree< value, boost::geometry::index::quadratic<16> > mTree;

    };
}
#endif
