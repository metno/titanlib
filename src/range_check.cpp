#include <vector>
#include <math.h>
#include "titanlib.h"
#include <assert.h>
#include <iostream>
#include <exception>
#include <boost/date_time/gregorian/gregorian_types.hpp>

using namespace titanlib;

double mean_temp(float lat, int unixtime);
double interpolate(int i1, int i2, float lat, float latitudes[], double numbers[]);
std::pair<int,int> find_between(float lat, float latitudes[], int len);


ivec titanlib::range_check(const vec& values,
        const vec& min,
        const vec& max) {

    // loop over all the lats/lons/elevs + value 
    // either min/max has length 1 or is the same length as the other vecs
    const int s = values.size();
    // assert that the max and min are either both size 1 or the size of values 
    if( (min.size() != s && min.size() != 1) || (max.size() != s && max.size() != 1) ) {
        throw std::runtime_error("Dimension mismatch");
    }

    ivec flags(s, 0);

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

    return flags;

}

ivec titanlib::range_check_climatology(const Points& points,
        const vec& values,
        int unixtime,
        const vec& pos,
        const vec& neg) {

    // loop over all the lats/lons/elevs + value 
    // either min/max has length 1 or is the same length as the other vecs
    const int s = points.size();
    const vec& lats = points.get_lats();
    const vec& lons = points.get_lons();
    const vec& elevs = points.get_elevs();

    if( lons.size() != s || elevs.size() != s || values.size() != s ) {
        throw std::runtime_error("Dimension mismatch");
    }
    if( (pos.size() != s && pos.size() != 1) || (neg.size() != s && neg.size() != 1) ) {
        throw std::runtime_error("Dimension mismatch");
    }

    ivec flags(s, 0);

    for(int i = 0; i < s; i++) {
        // leave the index to 0 if its the same max/min applied to everything
        // else same as loop
        int pos_i = (pos.size() == s) ? i : 0;
        int neg_i = (neg.size() == s) ? i : 0;

        // get best guess mean temp
        double t = mean_temp(lats[i], unixtime);
        std::cout << "mean t: " << t << "\n";
        double mean_pos = t + pos[pos_i];
        double mean_neg = t + neg[neg_i];

        // loop over the vectors and set the flags (0 = ok and 1 = bad)
        if(values[i] < mean_neg || values[i] > mean_pos) {
            flags[i] = 1;
        }
    }
    return flags;
}

double mean_temp(float lat, int unixtime) {

    // do some math to check if the value makes sense based on the lat/long
    //http://www-das.uwyo.edu/~geerts/cwx/notes/chap16/geo_clim.html        
    float mean_latitudes[] = {90,60,50,45,40,30,15,0,-15,-30,-35,-40,-45,-60,-90};
    double mean_temp[] = {-15,5,10,15,20,25,30,30,25,21,20,15,10,0,-25};

    // find what two points are closest?
    std::pair<int,int> p = find_between(lat, mean_latitudes, 15);
    //std::cout << " pair: " << latitudes[p.first] << " ";
    //std::cout << latitudes[p.second] << "\n";
    double mt = interpolate(p.first, p.second, lat, mean_latitudes, mean_temp); 

    // factor in season???
    //https://www.physics.byu.edu/faculty/christensen/physics%20137/Figures/Temperature/Effects%20of%20Latitude%20on%20Annual%20Temperature%20Range.htm
    float latitudes_var[] = {90,75,60,45,30,15,0,-15,-30,-45,-60,-75,90};
    double seasonal_variation[] = {40,32,30,23,13,3,0,4,7,6,11,26,31};

    // find what two points are closest?
    std::pair<int,int> p2 = find_between(lat, latitudes_var, 13);
    double sv = interpolate(p2.first, p2.second, lat, latitudes_var, seasonal_variation);

    boost::gregorian::date epoch(1970,1,1);
    boost::gregorian::date_duration diff(unixtime/86400);
    boost::gregorian::date newDate = epoch + diff;
    int date = newDate.year() * 10000 + newDate.month() * 100 + newDate.day();
    //std::cout << date << "\n";
    int m = newDate.month();
    //std::cout << "month: " << m << "\n";
   
    // given this month, adjust the mt based on sv
    if(m == 1) { // coldest
        return mt + (sv/2)*-1;
    }
    if(m == 12 || m == 2) { // second coldest
        return mt + (sv/2)*-0.6;
    }
    if(m == 11 || m == 3) { 
        return mt + (sv/2)*-0.3;
    }
    //if(m == 10 || m == 4) { // average??? 
    //    return mt + (sv/2)*0.6;
    //}
    if(m == 5 || m == 9) { 
        return mt + (sv/2)*0.3;
    }
    if(m == 6 || m == 8) { // second warmest
        return mt + (sv/2)*0.6;
    }
    if(m == 7) { // warmest
        return mt + (sv/2)*1;
    }
    return mt;
}
double interpolate(int i1, int i2, float lat, float latitudes[], double numbers[]) {
    
    // calculate percentage between the 2 and get mean temp based on this
    float diff = latitudes[i1] - lat;
    float space = latitudes[i1] - latitudes[i2];
    float percent = diff/space;
    //std::cout << " lat: " << lat;
    //std::cout << " percent: " << percent;
    //std::cout << "\n";

    float temp_diff = numbers[i2] - numbers[i1];
    //std::cout << " diff: ";
    //std::cout << temp_diff;
    //std::cout << "\n";
    double mt = numbers[i1] + (temp_diff*percent);
    //std::cout << "interpolated num: " << mt << "\n";
    return mt;
}
std::pair<int,int> find_between(float lat, float latitudes[], int len) {

    std::pair<float,float> p = std::make_pair(0,0);
    assert(lat <= 90);
    assert(lat >= -90);

    for(int i = 0; i < len; i++) {
        if(lat > latitudes[i]) {
            // then that means we are between this and the previous
            p = std::make_pair(i-1,i);
            break;
        }
    }
    return p;
}
