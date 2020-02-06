#include <vector>
#include <math.h>
#include "titanlib.h"
#include <assert.h>
#include <iostream>

double average_temp(float lat);
double find_closest(float lat, float latitudes[], int len);


bool titanlib::range_check(const fvec values,
        const fvec min,
        const fvec max,
        ivec& flags) {

    // loop over all the lats/lons/elevs + value 
    // either min/max has length 1 or is the same length as the other vecs
    const int s = values.size();
    // assert that the max and min are either both size 1 or the size of values 
    if( (min.size() != s && min.size() != 1) || (max.size() != s && max.size() != 1) ) { return false; }

    // resize the flags
    flags.resize(s, 0);

    for(int i = 0; i < s; i++) {
        // leave the index to 0 if its the same max/min applied to everything
        // else same as loop
        int min_i = (min.size() == s) ? i : 0;
        int max_i = (max.size() == s) ? i : 0;

        // loop over the vectors and set the flags (0 = ok and 1 = bad)
        if(values[i] < min[min_i] || values[i] > max[max_i]) {
            flags[i] = 1;
        }
    }

    return true;

}

bool range_check_climatology(const fvec lats,
        const fvec lons,
        const fvec elevs,
        const fvec values,
        const fvec min,
        const fvec max,
        ivec& flags) {

    // loop over all the lats/lons/elevs + value 
    // either min/max has length 1 or is the same length as the other vecs
    const int s = lats.size();
    if( lons.size() != s || elevs.size() != s || values.size() != s ) { return false; }
    if( (min.size() != s && min.size() != 1) || (max.size() != s && max.size() != 1) ) { return false; }

    flags.resize(s, 0);

    for(int i = 0; i < s; i++) {
        // leave the index to 0 if its the same max/min applied to everything
        // else same as loop
        int min_i = (min.size() == s) ? i : 0;
        int max_i = (max.size() == s) ? i : 0;

        // loop over the vectors and set the flags (0 = ok and 1 = bad)
        if(values[i] < min[min_i] || values[i] > max[max_i]) {
            flags[i] = 1;
        }
    
        // get best guess average temp
        average_temp(lats[i]);

    }

    return true;

}

double average_temp(float lat) {

    // do some math to check if the value makes sense based on the lat/long (and season?)
    //https://www.physics.byu.edu/faculty/christensen/physics%20137/Figures/Temperature/Effects%20of%20Latitude%20on%20Annual%20Temperature%20Range.htm
    //http://www-das.uwyo.edu/~geerts/cwx/notes/chap16/geo_clim.html        
    float latitudes[] = {90,60,50,45,40,30,15,0,-15,-30,-35,-40,-45,-60,-90};
    double mean_temp[] = {-15,5,10,15,20,25,30,30,25,21,20,15,10,0,-25};

    // find what two points are closest?
    find_closest(lat, latitudes, 15);

    return 1.0;
}
double find_closest(float lat, float latitudes[], int len) {
    
    for(int i = 0; i < len; i++) {
        if(lat > latitudes[i]) {
            // then that means we have gone past the closest one
            float diff1 = abs(lat - latitudes[i]);
            float diff2 = abs(lat - latitudes[i-1]);
            std::cout << "diff1: ";
            std::cout << diff1;
            std::cout << "diff2: ";
            std::cout << diff2;
        }
    }

}
