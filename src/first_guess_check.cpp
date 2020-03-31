#include <vector>
#include <iostream>
#include <math.h>
#include "titanlib.h"


int titanlib::first_guess_check(
    const fvec& values,
    const fvec& fg_values,
    const fvec& fg_neg,
    const fvec& fg_pos,
    ivec& flags) {

        const int s = values.size();
        //if( (lats.size() != s && lats.size() != 1) || (lons.size() != s && lons.size() != 1) ) { return false; }
        //f( (elevs.size() != s && elevs.size() != 1) || (values.size() != s && values.size() != 1) ) { return false; }
        if( (fg_values.size() != s && fg_values.size() != 1) ) { return 1; }
        if( (fg_neg.size() != s && fg_neg.size() != 1) || (fg_pos.size() != s && fg_pos.size() != 1) ) { return 1; }

        fvec min(s,0);
        fvec max(s,0);
        for(int i = 0; i < s; i++) { 
             int fgn_i = (fg_neg.size() == s) ? i : 0;
             int fgp_i = (fg_pos.size() == s) ? i : 0;
             min[i] = fg_values[i] - fg_neg[fgn_i];    
             max[i] = fg_values[i] + fg_pos[fgp_i]; 
        }
        
        // range check will resize the flags 
        return titanlib::range_check(values, min, max, flags);  
}
