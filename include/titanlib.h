#ifndef TITANLIB_H
#define TITANLIB_H
#include <iostream>
#include <vector>
#include <assert.h>
#include <libalglib/interpolation.h>
typedef std::vector<float> fvec;
typedef std::vector<double> dvec;
typedef std::vector<int> ivec;

/** Titanlib
*/
namespace titanlib {
    /** Spatial Consistency Test
      * @param lats vector of latitudes
      * @param flags output vector of flags
      */
    bool sct(const fvec lats,
            const fvec lons,
            const fvec elevs,
            const fvec values,
            int nmin,
            int nmax,
            int nminprof,
            float dzmin,
            float dhmin,
            float dz,
            const fvec t2pos,
            const fvec t2neg,
            const fvec eps2,
            fvec& sct,
            ivec& flags);

    bool range_check(const fvec values,
            const fvec min,
            const fvec max,
            ivec& flags);

    bool range_check_climatology(const fvec lats,
            const fvec lons,
            const fvec elevs,
            const fvec values,
            int unixtime,
            const fvec plus,
            const fvec minus,
            ivec& flags);

    bool buddy_check(const fvec lats,
            const fvec lons,
            const fvec elevs,
            const fvec values,
            const fvec distance_lim,
            const ivec priorities,
            const ivec buddies_min,
            const fvec thresholds,
            float diff_elev_max,
            bool adjust_for_elev_diff,
            ivec& flags,
            const ivec obs_to_check = ivec());

    /** Isolation check. Checks that a station is not located alone by itself
     *  @param lats vector of latitudes [deg]
     *  @param lons vector of longitudes [deg]
     *  @param nmin required number of observaions
     *  @param radius search radius [m]
     *  @param flags vector of return flags
     * */
    bool isolation_check(const fvec lats,
            const fvec lons,
            int nmin,
            float radius,
            ivec& flags);

    /** Isolation check with elevation. Checks that a station is not located alone by itself
     *  @param lats vector of latitudes [deg]
     *  @param lons vector of longitudes [deg]
     *  @param elevs vector of elevations [m]
     *  @param nmin required number of observaions
     *  @param radius search radius [m]
     *  @param dz vertical search radius [m]
     *  @param flags vector of return flags
     * */
    bool isolation_check(const fvec lats,
            const fvec lons,
            const fvec elevs,
            int nmin,
            float radius,
            float dz,
            ivec& flags);

    namespace util {
    /** Convert lat/lons to 3D cartesian coordinates with the centre of the earth as the origin
     *  @param lats vector of latitudes [deg]
     *  @param lons vector of longitudes [deg]
     *  @param x_coords vector of x-coordinates [m]
     *  @param y_coords vector of y-coordinates [m]
     *  @param z_coords vector of z-coordinates [m]
     * */
    bool convert_coordinates(const fvec& lats, const fvec& lons, fvec& x_coords, fvec& y_coords, fvec& z_coords);

    /** Same as above, but convert a single lat/lon to 3D cartesian coordinates
     *  @param lat latitude [deg]
     *  @param lon longitude [deg]
     *  @param x_coord x-coordinate [m]
     *  @param y_coord y-coordinate [m]
     *  @param z_coord z-coordinate [m]
     * */
    bool convert_coordinates(float lat, float lon, float& x_coord, float& y_coord, float& z_coord);
    }

    // ivec nearest_neighbours(const fvec& lats, const fvec& lons, float radius, float lat, float lon);

    // bool prioritize(const fvec values, const ivec priority, float distance, ivec& flags);
    class Dataset {
        public:
            Dataset(fvec ilats, fvec ilons, fvec ielevs, fvec ivalues);
            /** Perform the range check on the dataset
            * @param indices Only perform the test on these indices
            **/
            bool range_check(const fvec min, const fvec max, const ivec indices=ivec());
            bool range_check_climatology(int unixtime, const fvec plus, const fvec minus, const ivec indices=ivec());
            bool sct(int nmin, int nmax, int nminprof, float dzmin, float dhmin, float dz, const fvec t2pos, const fvec t2neg, const fvec eps2, fvec& sct, const ivec indices=ivec());
            bool buddy_check(const fvec distance_lim, const ivec priorities, const ivec buddies_min, const fvec thresholds, float diff_elev_max, bool adjust_for_elev_diff, const ivec obs_to_check, const ivec indices=ivec());

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
             * */
            int get_nearest_neighbour(float lat, float lon);

            /** Find all points with a radius
             *  @param lat Latitude of lookup-point
             *  @param lon Longitude of lookup-point
             *  @param radius Lookup radius [m]
             * */
            ivec get_neighbours(float lat, float lon, float radius);

            /** Find all points with a radius
             *  @param lat Latitude of lookup-point
             *  @param lon Longitude of lookup-point
             *  @param radius Lookup radius [m]
             *  @param distances Vector to store separation distances [m]
             * */
            ivec get_neighbours_with_distance(float lat, float lon, float radius, fvec& distances);

            /** Find the number of points within a radius
             *  @param lat Latitude of lookup-point
             *  @param lon Longitude of lookup-point
             *  @param radius Lookup radius [m]
             * */
            int get_num_neighbours(float lat, float lon, float radius);

            /** Find a set of nearest points
             *  @param lat Latitude of lookup-point
             *  @param lon Longitude of lookup-point
             *  @param num Number of points to find
             * */
            ivec get_closest_neighbours(float lat, float lon, int num);
        private:
            alglib::kdtree mTree;
            static alglib::real_1d_array ll2ar(float lat, float lon);

    };
}
#endif
