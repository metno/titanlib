#ifndef TITANLIB_H
#define TITANLIB_H
#include <iostream>
#include <vector>
#include <assert.h>
#include <boost/geometry.hpp>
#include <boost/geometry/geometries/point.hpp>
#include <boost/geometry/geometries/box.hpp>
#include <boost/geometry/index/rtree.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#ifdef _OPENMP
    #include <omp.h>
#endif

#define TITANLIB_VERSION "0.3.4.dev3"
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

    enum BackgroundType {
        VerticalProfile = 0,
        VerticalProfileTheilSen = 1,
        MeanOuterCircle = 2,
        MedianOuterCircle = 3,
        External = 4,
    };
 
    enum ConditionType {
        Eq = 0,
        Gt = 1,
        Geq = 2,
        Lt = 3,
        Leq = 4,
    };

    class Points;

    /** ****************************************
     * @name Spatial checks
     * Checks that consider the spatial properties of observations
     * *****************************************/ /**@{*/
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
     *  @param obs_to_check Observations that will be checked (since can pass in observations that will not be checked). 1=check the corresponding observation
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
            vec& rep,
            const ivec& obs_to_check=ivec());
    /** ****************************************
     * @name Spatial Checks
     * Validates observations based on spatial properties and first-guess values.
     * *****************************************/ /**@{*/
    /** Spatial Consistency Test with First Guess
     * @param points            Longitude, latitude, and elevation of observation locations
     * @param values            Observed values
     * @param background_values First-guess values at observation locations
     * @param values_min        Minimum acceptable observed value (set equal to values_max to ignore)
     * @param values_max        Maximum acceptable observed value (set equal to values_min to ignore)
     * @param num_min           Minimum required observations within the outer radius (must be > 1)
     * @param num_max           Maximum observations used for the test (must be > num_min)
     * @param inner_radius      Radius for flagging [m]
     * @param outer_radius      Radius for computing OI [m]
     * @param num_iterations    Maximum iterations (stops if no new flags are set)
     * @param min_horizontal_scale Minimum horizontal decorrelation length [m]
     * @param vertical_scale    Vertical decorrelation length [m]
     * @param pos               Allowed positive deviation
     * @param neg               Allowed negative deviation
     * @param eps2              Observation-to-background error variance ratio (e.g., 0.5 means observations are trusted twice as much as the background)
     * @param min_obs_var       Minimum observation error variance (reflects estimated representativeness error or expected observation uncertainty at min_horizontal_scale)
     * @param gross_error_score Gross error score per observation (higher values indicate a greater likelihood of measurement or large representativeness error)
     * @param obs_to_check      Observations to be checked (1 = check, 0 = ignore)
     * @return Flags indicating suspect observations
     */

    ivec sct_with_fg(const Points& points,
            const vec& values,
            const vec& background_values,
            float values_min,
            float values_max,
            int num_min,
            int num_max,
            float inner_radius,
            float outer_radius,
            int num_iterations,
            float min_horizontal_scale,
            float vertical_scale,
            const vec& pos,
            const vec& neg,
            const vec& eps2,
            const vec& min_obs_var,
            vec& gross_error_score,
            const ivec& obs_to_check=ivec());

    /** Spatial Consistency Test (SCT) - resistant to outliers
     *  @param points Input points
     *  @param values observed values to check (and/or to use)
     *  @param obs_to_check Observations that will be checked (since can pass in observations that will not be checked). 1=check the corresponding observation
     *  @param background_values external background value (not used if background_elab_type!=external)
     *  @param background_elab_type one of: VerticalProfile, VerticalProfileTheilSen, MeanOuterCircle, MedianOuterCircle, External
     *  @param num_min_outer Minimum number of observations inside the outer circle to compute SCT
     *  @param num_max_outer Maximum number of observations inside the outer circle used
     *  @param num_min_prof Minimum number of observations to compute vertical profile
     *  @param inner_radius Radius for flagging [m]
     *  @param outer_radius Radius for computing OI and background [m]
     *  @param num_iterations Number of SCT iterations
     *  @param min_elev_diff Minimum elevation difference to compute vertical profile [m]
     *  @param min_horizontal_scale Minimum horizontal decorrelation length [m]
     *  @param max_horizontal_scale Maximum horizontal decorrelation length [m]
     *  @param kth_closest_obs_horizontal_scale Number of closest observations to consider in the adaptive estimation of the horizontal decorrelation length
     *  @param vertical_scale Vertical decorrelation length [m]
     *  @param value_mina Minimum admissible value
     *  @param value_maxa Maximum admissible value
     *  @param value_minv Minimum valid value
     *  @param value_maxv Maximum valid value
     *  @param eps2 ratio between observation and background error variances
     *  @param tpos SCT-score threshold. Positive deviation allowed
     *  @param tneg SCT-score threshold. Negative deviation allowed
     *  @param debug Verbose output
     *  @param basic Which SCT-score should be used? Basic or more advanced, which takes into account the local variability. Boolean: true=basic; false=advanced
     *  @param scores SCT-score. The higher the score, the more likely is the presence of a gross measurement error
     *  @return flags
     */
    ivec sct_resistant( const Points& points,
                        const vec& values,
                        const ivec& obs_to_check,
                        const vec& background_values,
                        BackgroundType background_elab_type,
                        int num_min_outer,
                        int num_max_outer,
                        float inner_radius,
                        float outer_radius,
                        int num_iterations,
                        int num_min_prof,
                        float min_elev_diff,
                        float min_horizontal_scale,
                        float max_horizontal_scale,
                        int kth_closest_obs_horizontal_scale,
                        float vertical_scale,
                        const vec& value_mina,
                        const vec& value_maxa,
                        const vec& value_minv,
                        const vec& value_maxv,
                        const vec& eps2,
                        const vec& tpos,
                        const vec& tneg,
                        bool debug,
                        bool basic,
                        vec& scores);

    /**  Spatial Consistency Test for dichotomous (yes/no) variables
     *  @param points Input points
     *  @param values observed values to check (and/or to use)
     *  @param obs_to_check Observations that will be checked (since can pass in observations that will not be checked). 1=check the corresponding observation
     *  @param event_thresholds event thresholds (not used if background_elab_type!=external)
     *  @param condition one of: Eq, Gt, Geq, Lt, Leq
     *  @param num_min_outer Minimum number of observations inside the outer circle to compute SCT
     *  @param num_max_outer Maximum number of observations inside the outer circle used
     *  @param inner_radius Radius for flagging [m]
     *  @param outer_radius Radius for computing OI and background [m]
     *  @param num_iterations Number of SCT iterations
     *  @param min_horizontal_scale Minimum horizontal decorrelation length [m]
     *  @param max_horizontal_scale Maximum horizontal decorrelation length [m]
     *  @param kth_closest_obs_horizontal_scale Number of closest observations to consider in the adaptive estimation of the horizontal decorrelation length
     *  @param vertical_scale Vertical decorrelation length [m]
     *  @param test_thresholds observation-dependent trhesholds used in the SCT dual tests. The threshold is for the relative information content and it should tipically be a number between 0-1
     *  @param debug Verbose output
     *  @return flags
     */
    ivec sct_dual( const Points& points,
                   const vec& values,
                   const ivec& obs_to_check,
                   const vec& event_thresholds,
                   ConditionType condition,
                   int num_min_outer,
                   int num_max_outer,
                   float inner_radius,
                   float outer_radius,
                   int num_iterations,
                   float min_horizontal_scale,
                   float max_horizontal_scale,
                   int kth_closest_obs_horizontal_scale,
                   float vertical_scale,
                   const vec& test_thresholds,
                   bool debug);

     /** First Guess Test (FGT) - simplified (without OI) SCT
     *  @param points Input points
     *  @param values observed values to check (and/or to use)
     *  @param obs_to_check Observations that will be checked (since can pass in observations that will not be checked). 1=check the corresponding observation
     *  @param background_values external background value (not used if background_elab_type!=external)
     *  @param background_uncertainties uncertainty of the external background value (not used if background_elab_type!=external, optional when  background_elab_type=external)
     *  @param background_elab_type one of: vertical_profile, vertical_profile_Theil_Sen, mean_outer_circle, external
     *  @param num_min_outer Minimum number of observations inside the outer circle to compute FGT
     *  @param num_max_outer Maximum number of observations inside the outer circle used
     *  @param num_min_prof Minimum number of observations to compute vertical profile
     *  @param inner_radius Radius for flagging [m]
     *  @param outer_radius Radius for computing OI and background [m]
     *  @param num_iterations Number of FGT iterations
     *  @param min_elev_diff Minimum elevation difference to compute vertical profile [m]
     *  @param value_mina Minimum admissible value
     *  @param value_maxa Maximum admissible value
     *  @param value_minv Minimum valid value
     *  @param value_maxv Maximum valid value
     *  @param tpos FGT-score threshold. Positive deviation allowed
     *  @param tneg FGT-score threshold. Negative deviation allowed
     *  @param debug Verbose output
     *  @param basic Which SCT-score should be used? Basic or more advanced, which takes into account the local variability. Boolean: true=basic; false=advanced
     *  @param scores FGT-score. The higher the score, the more likely is the presence of a gross measurement error
     *  @return flags
     */
    ivec fgt( const Points& points,
              const vec& values,
              const ivec& obs_to_check,
              const vec& background_values,
              const vec& background_uncertainties,
              BackgroundType background_elab_type,
              int num_min_outer,
              int num_max_outer,
              float inner_radius,
              float outer_radius,
              int num_iterations,
              int num_min_prof,
              float min_elev_diff,
              const vec& value_mina,
              const vec& value_maxa,
              const vec& value_minv,
              const vec& value_maxv,
              const vec& tpos,
              const vec& tneg,
              bool debug,
              bool basic,
              vec& scores);

    /** Range check. Checks observation is within the ranges given
     *  @param values vector of observations
     *  @param min min allowed value
     *  @param max max allowed value
     *  @param flags vector of return flags
     */
    ivec range_check(const vec& values,
            const vec& min,
            const vec& max);

    ivec range_check_climatology(const Points& points,
            const vec& values,
            int unixtime,
            const vec& pos,
            const vec& neg);

    /** Buddy check. Compares a station to all its neighbours within a certain distance
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
    ivec buddy_check(const Points& points,
            const vec& values,
            const vec& radius,
            const ivec& num_min,
            float threshold,
            float max_elev_diff,
            float elev_gradient,
            float min_std,
            int num_iterations,
            const ivec& obs_to_check = ivec());

    ivec buddy_event_check(const Points& points,
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
     *  @param num_min required number of observations
     *  @param radius search radius [m]
     *  @param vertical_radius Vertical search radius [m]
     *  @param flags vector of return flags
     */
    ivec isolation_check(const Points& points,
            int num_min,
            float radius,
            float vertical_radius=MV);

    /** Isolation check, vectorized version.
     *  @param num_min required number of observations
     *  @param radius search radius [m]
     *  @param vertical_radius Vertical search radius [m]. Use empty array to disable.
     *  @param flags vector of return flags
     */
    ivec isolation_check(const Points& points,
            const ivec& num_min,
            const vec& radius,
            const vec& vertical_radius=vec());

    /** Duplicate check. Checks Flags duplicate stations. Keeps the first one when duplicates.
     *  @param radius Stations within this radius are considered duplicates [m]
     *  @param vertical_range Stations must also be within this vertical range to be considered duplicates [m]. Disable if MV.
     */
    ivec duplicate_check(const Points& points, float radius, float vertical_range=titanlib::MV);

    ivec metadata_check(const Points& points, bool check_lat=true, bool check_lon=true, bool check_elev=true, bool check_laf=true);

    /** ****************************************
     * @name Timeserie methods
     * Functions that operate on timeseries of observations
     * *****************************************/ /**@{*/
    /** Method by McCarthy 1973
      * https://doi.org/10.1175/1520-0450(1973)012%3C0211:AMFCAT%3E2.0.CO;2
    */
    vec lag_reduction_filter(const vec& times, const vec& values, float a=0.5, float b=0.5, float k1=0.25, float k2=0.25, int n=10);

    /** ****************************************
     * @name OpenMP settings
     * Functions that configure OpenMP
     * *****************************************/ /**@{*/
    /** Set the number of OpenMP threads to use. Overwrides OMP_NUM_THREAD env variable. */
    void set_omp_threads(int num);

    /** Get the number of OpenMP threads currently set */
    int get_omp_threads();

    /** Sets the number of OpenMP threads to 1 if OMP_NUM_THREADS undefined */
    void initialize_omp();

    /** ****************************************
     * @name Utilities
     * Helper functions
     * *****************************************/ /**@{*/

    /**
     * @return Titanlib version
     */
    std::string version();

    /**
     * @return The current UTC time (in seconds since 1970-01-01 00:00:00Z)
     */
    double clock();
    bool is_valid(float value);

    /** Compute a vertical profile based on input data
     *  @param elevs vector of input elevations [m]
     *  @param oelevs vector of output elevatios [m]
     *  @param values observation vector [degrees C]
     *  @param num_min_prof minimum number of observations needed to do a full profile
     *  @param min_elev_diff use no elevation profile if elevation difference among stations is less than this [m]
     *  @param debug If true, print debug messages
     */
    vec compute_vertical_profile(const vec& elevs, const vec& oelevs, const vec& values, int num_min_prof, double min_elev_diff, bool debug=false);

    vec compute_vertical_profile_Theil_Sen(const vec& elevs, const vec& oelevs, const vec& values, int num_min_prof, double min_elev_diff, bool debug);

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
    float find_k_closest(const vec& array, int k);
    template <class T> T subset(const T& array, const ivec& indices) {
        if(array.size() == 0)
            return array;
        if(indices.size() == 0) {
            return array;
        }

        T new_array(indices.size());
        for(int i = 0; i < indices.size(); i++) {
            new_array[i] = array[indices[i]];
        }
        return new_array;
    }
    Points subset(const Points& input, const ivec& indices);
    template <class T> void unsubset(const T& array, T& orig_array, const ivec& indices) {
        orig_array.clear();
        orig_array.resize(indices.size());
        assert(array.size() == indices.size());
        for(int i = 0; i < indices.size(); i++) {
            orig_array[indices[i]] = array[i];
        }
    }


    /** Background (first guess) calculations at observation locations
    *  @param background_type background elaboration type
    *  @param elevs elevations (m above mean sea level)
    *  @param values observed values
    *  @param num_min_prof minimum number of observations needed to compute a vertical profile different from the pseudo-adiabatic lapse rate (-0.0065 degC/m)
    *  @param min_elev_diff when computing vertical profiles, the elevation difference between the 95th percentile and 5th percentile must be greater than this parameter. Otherwise use the pseudo-adiabatic lapse rate
    *  @param value_minp minimum plausible value
    *  @param value_maxp maximum plausible value
    *  @param external_background_values external background values (used when background_elab_type=external)
    *  @param indices_global_outer vector of positions of matches of (vector with the observations belonging to the outer circle) in the (global vector) (dimension is the number of observations in the outer circle)
    *  @param debug verbose mode (true / false)
    */
    vec background(const vec& elevs, const vec& values, int num_min_prof, float min_elev_diff, float value_minp, float value_maxp, BackgroundType background_type, const vec& external_background_values, const ivec& indices_global_outer, bool debug);

    bool invert_matrix (const boost::numeric::ublas::matrix<float>& input, boost::numeric::ublas::matrix<float>& inverse);

    bool set_indices( const ivec& indices_global_outer_guess, const ivec& obs_test, const ivec& dqcflags, const vec& dist_outer_guess, float inner_radius, int test_just_this, ivec& indices_global_outer, ivec& indices_global_test, ivec& indices_outer_inner, ivec& indices_outer_test, ivec& indices_inner_test);
    template<class T1, class T2> struct sort_pair_first {
        bool operator()(const std::pair<T1,T2>&left, const std::pair<T1,T2>&right) {
            return left.first < right.first;
        };
    };
    /**@}*/

    /** ****************************************
     * @name SWIG testing functions
     * Functions for testing the SWIG interface. Not useful for any other purpose.
     * *****************************************/ /**@{*/
    /** Required for SWIG only */
    float* test_array(float* v, int n);

    void test_not_implemented_exception();
    /**@}*/
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
            KDTree(CoordinateType type=Geodetic) {mType=type;};

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
             *  @param indices Only perform the test on these indices. Use ivec(1, -1) to indicate all stations
             */
            void range_check(const vec& min, const vec& max, const ivec& indices=ivec(1, -1));
            void range_check_climatology(int unixtime, const vec& pos, const vec& neg, const ivec& indices=ivec(1, -1));
            void sct(int num_min, int num_max, 
                    float inner_radius, float outer_radius, 
                    int num_iterations, int num_min_prof, 
                    float min_elev_diff, float min_horizontal_scale, float vertical_scale, 
                    const vec& t2pos, const vec& t2neg, const vec& eps2, vec& prob_gross_error, 
                    vec& rep, const ivec& indices=ivec(1, -1));
            void sct_dual(const vec& event_thresholds, 
                                ConditionType condition, 
                                int num_min, int num_max,
                                float inner_radius, float outer_radius,
                                int num_iterations,
                                float min_horizontal_scale, float max_horizontal_scale,
                                int kth_closest_obs_horizontal_scale,
                                float vertical_scale,
                                const vec& test_thresholds,
                                bool debug,
                                const ivec& obs_to_check=ivec(), const ivec& indices=ivec(1, -1));
            void buddy_check(const vec& radius, const ivec& num_min, float threshold, float max_elev_diff, float elev_gradient, float min_std, int num_iterations,
                             const ivec& obs_to_check=ivec(), const ivec& indices=ivec(1, -1));
            void buddy_event_check(const vec& radius, const ivec& num_min, float event_threshold, float threshold, float max_elev_diff, float elev_gradient, int num_iterations,
                             const ivec& obs_to_check=ivec(), const ivec& indices=ivec(1, -1));
            void isolation_check(int num_min, float radius, float vertical_radius=MV, const ivec& indices=ivec(1, -1));
            void isolation_check(const ivec& num_min, const vec& radius, const vec& vertical_radius=vec(), const ivec& indices=ivec(1, -1));
            void duplicate_check(float radius, float vertical_range=titanlib::MV, const ivec& indices=ivec(1, -1));
            void dem_check(const vec& dem, float max_elev_diff);
            void external_check(const ivec& flags);
            void metadata_check(bool check_lat=true, bool check_lon=true, bool check_elev=true, bool check_laf=true, const ivec& indices=ivec(1, -1));

            Points get_points() const;
            vec get_values() const;
            ivec get_flags() const;
            void set_values(vec ivalues);
            void set_flags(ivec ivalues);
            void set_points(Points ipoints);

            Points points;
            vec values;
            ivec flags;
        private:
            /* Get a subset of array where flags=0 and for indices. This function is different than
             * titanlib::subset in that it takes into account flags.
             * @param array Vector of values
             * @param indices Vector of indices
             * @return Vector of values
            */
            template <class T> T subset_valid(const T& array, const ivec& indices=ivec(1, -1)) {
                if(array.size() != 1 && array.size() != flags.size()) {
                    std::stringstream ss;
                    ss << "Array (" << array.size() << ") must be either size 1 or the same as dataset size (" << flags.size() << ")";
                    throw std::invalid_argument(ss.str());
                }
                if(indices.size() == 0)
                    return T();

                ivec indices0 = indices;
                if(indices.size() == 1 && indices[0] == -1) {
                    // Use all indices
                    indices0.clear();
                    indices0.resize(flags.size());
                    for(int i = 0; i < flags.size(); i++)
                        indices0[i] = i;
                }
                T new_array;
                new_array.reserve(indices0.size());
                for(int i = 0; i < indices0.size(); i++) {
                    int index = indices0[i];
                    if(flags[index] == 0) {
                        if(array.size() == 1)
                            // Array is broadcast
                            new_array.push_back(array[0]);
                        else {
                            assert(index < array.size());
                            new_array.push_back(array[index]);
                        }
                    }
                }
                return new_array;
            }
            /*  Same as T subset_valid, except for Points */
            Points subset_valid(const Points& input, const ivec& indices=ivec(1, -1));

            ivec get_unflagged_indices();

            // Create a vector with results, but only for locations that have not been flagged
            // To broadcase the values to point, use a singleton array
            // An empty array will return an empty array
            template <class T> T get_unflagged(const T& array) {
                if(array.size() == 0)
                    return array;

                if(array.size() > 1 && array.size() < flags.size()) {
                    throw std::invalid_argument("Array is shorter than flags");
                }

                ivec indices = get_unflagged_indices();
                T output(indices.size());
                for(int i = 0; i < indices.size(); i++) {
                    if(array.size() > 1) {
                        int index = indices[i];
                        assert(index < array.size());
                        output[i] = array[index];
                    }
                    else
                        output[i] = array[0];
                }
                return output;
            }
            Points get_unflagged_points();

            /* Same as T subset_valid, except the subsetted indices are returned */
            ivec subset_valid(const ivec& indices);

            /* Merge the flags into the objects flags. If existing flags are 1 but new_flags are
             * 0, the existing flags are not updated.
             * @param new_flags Vector of new flags
             * @param indices The location indices that flags are valid for
            */
            void merge_simple(const ivec& new_flags, ivec indices);

            /* Merge new flags into the object's flags. If existing flags are 1 but new_flags are
             * 0, the existing flags are not updated.
             * @param new_flags Vector of new flags
             * @param indices The location indices that flags are valid for
            */
            void merge(const ivec& new_flags, ivec indices);
    };
    class not_implemented_exception: public std::logic_error
    {
        public:
            not_implemented_exception() : std::logic_error("Function not yet implemented") { };
    };

}
#endif
