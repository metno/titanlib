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

#define TITANLIB_VERSION "0.2.0b3"
#define __version__ TITANLIB_VERSION

namespace titanlib {
    /** **************************************
     * @name Short-hand notation for vectors of different dimensions sizes
     * ***************************************/ /**@{*/
    // Preferred vector types
    typedef std::vector<int> ivec;
    typedef std::vector<float> vec;
    typedef std::vector<double> dvec;
    typedef std::vector<vec> vec2;
    /**@}*/

    /** **************************************
     * @name Constants
     * Functions that assimilate observations onto a gridded background
     * ***************************************/ /**@{*/
    /** Missing value indicator */
    static const float MV = NAN;
    /** Missing value indicator in gridpp command-line tool */
    static const float MV_CML = -999;
    /** Mathematical constant pi */
    static const float pi = 3.14159265;
    /** Radius of the earth [m] */
    static const double radius_earth = 6.378137e6;
    /**@}*/

    /** Types of coordinates for position of points */
    enum CoordinateType {
        Geodetic = 0,      /**< Latitude and longitude */
        Cartesian = 1,     /**< X and Y */
    };

    class Points;

    /**
     * @return Titanlib version
     */
    std::string version();

    /** Spatial Consistency Test
     *  @param num_min_prof Minimum number of observations to compute vertical profile
     *  @param inner_radius Radius for flagging [m]
     *  @param outer_radius Radius for computing OI and background [m]
     *  @param min_elev_diff Minimum elevation difference to compute vertical profile [m]
     *  @param min_horizontal_scale Minimum horizontal decorrelation length [m]
     *  @param vertical_scale Vertical decorrelation length [m]
     *  @param pos Positive deviation allowed
     *  @param neg Negative deviation allowed
     *  @param eps2
     *  @param prob_gross_error Probability of gross error for each observation
     *  @param rep Coefficient of representativity
     *  @return flags
     */
    ivec sct(const Points& points,
            const vec& values,
            int num_min,
            int num_max,
            float inner_radius,
            float outer_radius,
            int num_iterations,
            int num_min_prof,
            float min_elev_diff,
            float min_horizontal_scale,
            float vertical_scale,
            const vec& pos,
            const vec& neg,
            const vec& eps2,
            vec& prob_gross_error,
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
     *  @param num_min the minimum number of buddies a station can have
     *  @param threshold the threshold for flagging a station
     *  @param max_elev_diff the maximum difference in elevation for a buddy (if negative will not check for heigh difference)
     *  @param elev_gradient linear elevation gradient with height. For temperature, use something like -0.0065.
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
            const ivec& num_min,
            float threshold,
            float max_elev_diff,
            float elev_gradient,
            float min_std,
            int num_iterations,
            const ivec& obs_to_check = ivec());

    ivec buddy_event_check(const vec& lats,
            const vec& lons,
            const vec& elevs,
            const vec& values,
            const vec& radius,
            const ivec& num_min,
            float event_threshold,
            float threshold,
            float max_elev_diff,
            float elev_gradient,
            int num_iterations,
            const ivec& obs_to_check = ivec());

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
     *  @param lats Vector of latitudes [deg]
     *  @param lons Vector of longitudes [deg]
     *  @param elevs Vector of elevations [m]
     *  @param num_min Required number of observations
     *  @param radius Search radius [m]
     *  @param vertical_radius Vertical search radius [m]
     *  @param flags Vector of return flags
     */
    ivec isolation_check(const vec& lats,
            const vec& lons,
            const vec& elevs,
            int num_min,
            float radius,
            float vertical_radius);

    /** Method by McCarthy 1973 */
    vec lag_reduction_filter(const vec& times, const vec& values, float a=0.5, float b=0.5, float k1=0.25, float k2=0.25, int n=10);

    /** ****************************************
     * @name OpenMP settings
     * Functions that configure OpenMP
     * *****************************************/ /**@{*/
    /** Set the number of OpenMP threads to use. Overwrides OMP_NUM_THREAD env variable. */
    void set_omp_threads(int num);

    /** Sets the number of OpenMP threads to 1 if OMP_NUM_THREADS undefined */
    void initialize_omp();

    /** ****************************************
     * @name Utilities
     * Helper functions
     * *****************************************/ /**@{*/
    /**
     * @return The current UTC time (in seconds since 1970-01-01 00:00:00Z)
     */
    double clock();
    bool is_valid(float value);

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

    vec interpolate_to_points(const vec2& input_lats, const vec2& input_lons, const vec2& input_values, const vec& output_lats, const vec& output_lons);
    float deg2rad(float deg);
    float calc_distance(float lat1, float lon1, float lat2, float lon2);
    float calc_distance(float x0, float y0, float z0, float x1, float y1, float z1);

    float compute_quantile(double quantile, const vec& array);
    vec subset(const vec& input, const ivec& indices);
    Points subset(const Points& input, const ivec& indices);
    /**@}*/

    /** ****************************************
     * @name SWIG testing functions
     * Functions for testing the SWIG interface. Not useful for any other purpose.
     * *****************************************/ /**@{*/
    /** Required for SWIG only */
    float* test_array(float* v, int n);
    /**@}*/

    void test_not_implemented_exception();
    // ivec nearest_neighbours(const vec& lats, const vec& lons, float radius, float lat, float lon);

    // bool prioritize(const vec& values, const ivec& priority, float distance, ivec& flags);

    class Point {
        public:
            /** Constructor
              * @param lat: Latitude coordinate
              * @param lon: Longitude coordinate
              * @param elev: Elevation
              * @param laf: Land area fraction (between 0 and 1)
              * @param type: Coordinate type for lat and lon
            */
            Point(float lat, float lon, float elev=MV, float laf=MV, CoordinateType type=Geodetic);
            float lat;
            float lon;
            float elev;
            float laf;
            CoordinateType type;
    };
    class KDTree {
        public:
            KDTree(vec lats, vec lons, CoordinateType type=Geodetic);
            KDTree& operator=(KDTree other);
            KDTree(const KDTree& other);
            KDTree() {};

            /** Find single nearest points
             *  @param lat Latitude of lookup-point
             *  @param lon Longitude of lookup-point
             * */
            int get_nearest_neighbour(float lat, float lon, bool include_match=true) const;

            /** Find all points with a radius
             *  @param lat Latitude of lookup-point
             *  @param lon Longitude of lookup-point
             *  @param radius Lookup radius [m]
             * */
            ivec get_neighbours(float lat, float lon, float radius, bool include_match=true) const;

            /** Find all points with a radius
             *  @param lat Latitude of lookup-point
             *  @param lon Longitude of lookup-point
             *  @param radius Lookup radius [m]
             *  @param distances Vector to store separation distances [m]
             * */
            ivec get_neighbours_with_distance(float lat, float lon, float radius, vec& distances, bool include_match=true) const;

            /** Find the number of points within a radius
             *  @param lat Latitude of lookup-point
             *  @param lon Longitude of lookup-point
             *  @param radius Lookup radius [m]
             * */
            int get_num_neighbours(float lat, float lon, float radius, bool include_match=true) const;

            /** Find a set of nearest points
             *  @param lat Latitude of lookup-point
             *  @param lon Longitude of lookup-point
             *  @param num Number of points to find
             * */
            ivec get_closest_neighbours(float lat, float lon, int num, bool include_match=true) const;


            /** Convert lat/lons to 3D cartesian coordinates with the centre of the earth as the origin
             *  @param lats vector of latitudes [deg]
             *  @param lons vector of longitudes [deg]
             *  @param x_coords vector of x-coordinates [m]
             *  @param y_coords vector of y-coordinates [m]
             *  @param z_coords vector of z-coordinates [m]
             * */
            bool convert_coordinates(const vec& lats, const vec& lons, vec& x_coords, vec& y_coords, vec& z_coords) const;

            /** Same as above, but convert a single lat/lon to 3D cartesian coordinates
             *  @param lat latitude [deg]
             *  @param lon longitude [deg]
             *  @param x_coord x-coordinate [m]
             *  @param y_coord y-coordinate [m]
             *  @param z_coord z-coordinate [m]
             * */
            bool convert_coordinates(float lat, float lon, float& x_coord, float& y_coord, float& z_coord) const;
            static float deg2rad(float deg);
            static float rad2deg(float deg);
            static float calc_distance(float lat1, float lon1, float lat2, float lon2, CoordinateType type=Geodetic);
            static float calc_distance(float x0, float y0, float z0, float x1, float y1, float z1);
            static float calc_distance_fast(float lat1, float lon1, float lat2, float lon2, CoordinateType type=Geodetic);
            static float calc_distance_fast(const Point& p1, const Point& p2);
            vec get_lats() const;
            vec get_lons() const;
            int size() const;
            CoordinateType get_coordinate_type() const;
        protected:
            typedef boost::geometry::model::point<float, 3, boost::geometry::cs::cartesian> point;
            typedef std::pair<point, unsigned> value;
            typedef boost::geometry::model::box<point> box;
            boost::geometry::index::rtree< value, boost::geometry::index::quadratic<16> > mTree;
            vec mLats;
            vec mLons;
            CoordinateType mType;

            struct within_radius {
                public:
                    within_radius(point p, float radius);
                    bool operator()(value const& v) const;
                private:
                    float radius;
                    point p;
            };
            struct is_not_equal {
                public:
                    is_not_equal(point p);
                    bool operator()(value const& v) const;
                private:
                    point p;
            };
    };
    class Points  {
        public:
            Points();
            /** Initialize a new grid
             *  @param lats: vector of latitudes [degrees]
             *  @param lons: vector of longitudes [degrees]
             *  @param elevs: vector of elevations [m]
             *  @param lafs: vector of land area fractions [1]
             *  @param type: Coordinate type
            */
            Points(vec lats, vec lons, vec elevs=vec(), vec lafs=vec(), CoordinateType type=Geodetic);
            Points(KDTree tree, vec elevs=vec(), vec lafs=vec());
            Points& operator=(Points other);
            Points(const Points& other);
            // Returns -1 if there are no neighbours
            int get_nearest_neighbour(float lat, float lon, bool include_match=true) const;
            ivec get_neighbours(float lat, float lon, float radius, bool include_match=true) const;
            ivec get_neighbours_with_distance(float lat, float lon, float radius, vec& distances, bool include_match=true) const;
            int get_num_neighbours(float lat, float lon, float radius, bool include_match=true) const;
            ivec get_closest_neighbours(float lat, float lon, int num, bool include_match=true) const;

            vec get_lats() const;
            vec get_lons() const;
            vec get_elevs() const;
            vec get_lafs() const;
            int size() const;
            CoordinateType get_coordinate_type() const;
        private:
            KDTree mTree;
            vec mLats;
            vec mLons;
            vec mElevs;
            vec mLafs;
    };

    /** Represents point and their observed values */
    class Dataset {
        public:
            Dataset(Points points, vec ivalues);
            /** Perform the range check on the dataset
             *  @param indices Only perform the test on these indices
             */
            void range_check(const vec& min, const vec& max, const ivec& indices=ivec());
            void range_check_climatology(int unixtime, const vec& pos, const vec& neg, const ivec& indices=ivec());
            void sct(int num_min, int num_max, float inner_radius, float outer_radius, int num_iterations, int num_min_prof, float min_elev_diff, float min_horizontal_scale, float vertical_scale, const vec& t2pos, const vec& t2neg, const vec& eps2, vec& sct, vec& rep, const ivec& indices=ivec());
            void buddy_check(const vec& radius, const ivec& num_min, float threshold, float max_elev_diff, float elev_gradient, float min_std, int num_iterations, const ivec& obs_to_check, const ivec& indices=ivec());
            void buddy_event_check(const vec& radius, const ivec& num_min, float event_threshold, float threshold, float max_elev_diff, float elev_gradient, int num_iterations, const ivec& obs_to_check = ivec(), const ivec& indices=ivec());
            void isolation_check(int num_min, float radius, float vertical_radius);

            vec lats;
            vec lons;
            vec elevs;
            Points points;
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
    class not_implemented_exception: public std::logic_error
    {
        public:
            not_implemented_exception() : std::logic_error("Function not yet implemented") { };
    };

}
#endif
