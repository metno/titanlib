#include <vector>
#include <iostream>
#include <math.h>
#include <exception>
#include "titanlib.h"


ivec titanlib::first_guess_check(
    const vec& values,
    const vec& first_guess,
    const vec& neg,
    const vec& pos) {

        const int s = values.size();
        //if( (lats.size() != s && lats.size() != 1) || (lons.size() != s && lons.size() != 1) ) { return false; }
        //f( (elevs.size() != s && elevs.size() != 1) || (values.size() != s && values.size() != 1) ) { return false; }
        if( (first_guess.size() != s && first_guess.size() != 1) ) {
            throw std::runtime_error("Dimension mismatch");
        }
        if( (neg.size() != s && neg.size() != 1) || (pos.size() != s && pos.size() != 1) ) {
            throw std::runtime_error("Dimension mismatch");
        }

        vec min(s,0);
        vec max(s,0);
        for(int i = 0; i < s; i++) { 
             int fgn_i = (neg.size() == s) ? i : 0;
             int fgp_i = (pos.size() == s) ? i : 0;
             min[i] = first_guess[i] - neg[fgn_i];    
             max[i] = first_guess[i] + pos[fgp_i]; 
        }
        
        // range check will resize the flags 
        return titanlib::range_check(values, min, max);  
}
