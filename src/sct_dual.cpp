#include <vector>
#include <math.h>
#include "titanlib.h"
#include <assert.h>
#include <iostream>
#include <numeric>
#include <exception>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/vector_proxy.hpp>
#include <boost/numeric/ublas/triangular.hpp>
#include <boost/numeric/ublas/lu.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/ublas/matrix.hpp>

using namespace titanlib;

// SCT_dual within the inner circle

bool sct_dual_core( const vec& lats, const vec& lons, const vec& elevs, const vec& w, const vec& t, float Dh_min, float Dh_max, int kth_close, float Dz, const vec& eps2, const ivec& indices_global_test, const ivec& indices_outer_inner, const ivec& indices_outer_test, const ivec& indices_inner_test, bool debug, float na, bool set_flag0, int& thrown_out, ivec& flags); 


//=============================================================================

//+ SCT_dual - Spatial Consistency Test for dichotomous (yes/no) variables
ivec titanlib::sct_dual( const Points& points,
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
                         bool debug) {
/*

 -+-+ Spatial Consistency Test for dichotomous (yes/no) variables +-+-

 + Description: 
    Flag observations that are (likely) affected by gross measurement errors (GE) based on neighbouring observations. The observations are transformed into dichotomous predictions of the type, "yes, an event happened", or "no, an event did not happen". The term "dual" is related to the fact that observations are divided into two classes. Observations are affected by GEs if their values are (a) [not related to the actual atmospheric state] OR (b) [affected by such large representativeness errors (REs) that they are difficult to reconstruct using the neighbouring observations].
  For simplicity, we will refer to observations affected by GE as bad observations, viceversa good observations are those not affected by GEs. 

 + Reference for some concepts, such as SCT and Integral Data Influence (IDI): 
    Lussana, C., Uboldi, F. and Salvati, M.R. (2010), A spatial consistency test for surface observations from mesoscale meteorological networks. Q.J.R. Meteorol. Soc., 136: 1075-1088. doi:10.1002/qj.622 (Reference Abbreviation: LUS10)
 See also the wiki-pages https://github.com/metno/titanlib/wiki/Spatial-consistency-test

 + Algorithm: 
    The big task of applying the SCT_dual over the whole observational network is splitted into several smaller tasks.
  The observations are divided into the two classes w=1) "yes, an event happened", or w=0) "no, an event did not happen"
  The SCT_dual is applied many times in sequence. 
  Each time, the SCT_dual is applied over a small region centered on an observation location (i.e. the centroid observation).
  Based on the centroid, we define an inner circle and an outer circle.
  The outer circle is used to provide the boundary conditions.
  The observations to test are those within the inner circle that have not been tested yet.
  The quality of those observations is determined by comparing their classes against the leave-one-out predicted classes considering all observations within the outer circle.
  See also the comments of the function "sct_dual_core" for more details.
  
  The application of the sequence of SCT_duals over the whole observational network is repeated several times, until no GEs are found within the observations tested. 
  In fact, the core of SCT_dual is to test a single observation against an independent estimate that is obtained assuming all the neighbouring observations are good ones. As a consequence, a bad observation can influence the results of some tests on neighbouring observations before being rejected.
  We apply special care to avoid "false positives" (good flagged as bad).
  For each SCT_dual over a centroid, we flag as bad only the observation showing the worst result.
  In addition, an extra-loop is done over the bad observations where only the good ones are considered. We may "save" bad observations.
  Analogously, to avoid "false negatives" (bad flagged as good), the last loop is done over the bad observations where only the good ones are considered. We may "save" bad observations during this last loop.

  (*) Notes on the data structures (see also comments in set_indices): we have defined a hierarchy of data structures
    global -> outer -> inner -> test
   global vectors are the input vectors. The other data structures refer to a specific iteration when a centroid observation has been identified. outer are the vectors of data on the outer circle, not yet flagged as bad. inner are the data on the inner circle, not yet flagged as bad. test data are the observations in the inner circle and still to test.
   We use the same counter variables for the data structures: g -> global; i,j -> outer; l -> inner; m -> test
   There are vectors of indices matching the different structures. 
   For example, indices_outer_inner matches the outer and the inner vectors. indices_outer_inner is a vector of the dimension of the number of observations in the inner circle. The instruction i=indices_outer_inner[l] maps the l-th element of inner-type vectors as the i-th element of outer-type vectors.

 (*) Notes on event_thresholds. Different thresholds can be used for different observations. The thresholds vectors have the same dimensions of the observations vectors (global).

 + Pseudo-code:
 
  Loop over the number of iterations prescribed
    Loop over the p observations (observation index is "curr"). 
      Decide if the observation is suitable for a centroid
      Define outer/inner cirles, observations to test and set all the vector of indices
      Check if the observation is isolated
      If all the observations in the outer circle belong to the same class, then further steps are not needed and the observations to test are flagged as good ones
      SCT_dual over the inner circle, either flag the worse observation as bad or flag all of them as good ones
      NOTE: observations can be flagged as good ones only after the first iteration (after rejecting the worst)
  
  Loop over the observations without QC flags (nor good neither bad) and use each of them as a centroid

  Loop over the bad observations and use each of them as a centroid
  
 + Returned values:

  flags. -999 = not checked; 0 = passed (good); 1 = failed (bad); 11 = isolated (<2 inside inner); 12 = isolated (<num_min_outer inside outer)

*/

    double s_time = titanlib::clock();

    const vec& lats = points.get_lats();
    const vec& lons = points.get_lons();
    const vec& elevs = points.get_elevs();
    const int p = values.size();
    if( lats.size() != p || lons.size() != p || elevs.size() != p || values.size() != p)
        throw std::runtime_error("Dimension mismatch");
    if(num_min_outer < 2)
        throw std::invalid_argument("num_min_outer must be > 1");
    if(num_max_outer < num_min_outer)
        throw std::invalid_argument("num_max_outer must be > num_min_outer");
    if(num_iterations < 1)
        throw std::invalid_argument("num_iterations must be >= 1");
    if(min_horizontal_scale <= 0)
        throw std::invalid_argument("min_horizontal_scale must be > 0");
    if(kth_closest_obs_horizontal_scale <= 0)
        throw std::invalid_argument("kth_closest_obs_horizontal_scale must be > 0");
    if(vertical_scale <= 0)
        throw std::invalid_argument("vertical_scale must be > 0");
    if(inner_radius < 0)
        throw std::invalid_argument("inner_radius must be >= 0");
    if(outer_radius < inner_radius)
        throw std::invalid_argument("outer_radius must be >= inner_radius");

    // initializations
    float na = -999.; // code for "Not Available". Any type of missing data
    int flag_not_tested_inner;
    int flag_not_tested_outer;
    ivec flags( p, na);
    ivec obs_test( p, 1);
    vec w( p, 0);
    vec eps2( p, 0.1);
    vec r( p, na);
    vec t( p, na);

    bool set_all_good = false;
    bool accept_not_tested = true;
    if(accept_not_tested){
        std::cout << "Set flags for values that could not be tested as good " << p << std::endl;
        flag_not_tested_inner = 0;
        flag_not_tested_outer = 0;
    }
    else{
        std::cout << "Set flags for values that could not be tested due to lack of points in the inner/outer radius to 11/12." << p << std::endl;
        flag_not_tested_inner = 11;
        flag_not_tested_outer = 12;
    }
    
    if (obs_to_check.size() == p) {
        for(int g=0; g<p; g++) 
          obs_test[g] = obs_to_check[g];
    }

    if (event_thresholds.size() == p) {
        for(int g=0; g<p; g++) 
          r[g] = event_thresholds[g];
    } else {
        for(int g=0; g<p; g++) 
          r[g] = event_thresholds[0];
    }

    if (test_thresholds.size() == p) {
        for(int g=0; g<p; g++) 
          t[g] = test_thresholds[g];
    } else {
        for(int g=0; g<p; g++) 
          t[g] = test_thresholds[0];
    }

    for(int g=0; g<p; g++) {
      if ( condition == titanlib::Eq) {
        if ( values[g] == r[g]) w[g] = 1;
      } else if ( condition == titanlib::Geq) {
        if ( values[g] >= r[g]) w[g] = 1;
      } else if ( condition == titanlib::Gt) {
        if ( values[g] >  r[g]) w[g] = 1;
      } else if ( condition == titanlib::Leq) {
        if ( values[g] <= r[g]) w[g] = 1;
      } else if ( condition == titanlib::Lt) {
        if ( values[g] <  r[g]) w[g] = 1;
      }
    }

    if(debug) std::cout << "==================DEBUG MODE ============================= " << std::endl;
    if(debug) std::cout << "Number of observations to test is " << p << std::endl;

    // KDtree has to do with fast computation of distances
    titanlib::KDTree tree(lats, lons);

    //
    if(debug) {
        std::cout << "g lats lons elevs obs:" << std::endl;
        for(int g=0; g < p; g++) {
            std::cout << g << " " << lats[g] << " " << lons[g] << " " << elevs[g] << " " << values[g] << std::endl;
        }
    }

    // SCT_dual iterations
    for(int iteration = 0; iteration < num_iterations; iteration++) {

        if(debug) std::cout << " +++++ Iteration " << iteration << " ++++++++++++++++" << std::endl;

        double s_time0 = titanlib::clock();
        
        // reset this number each loop (this is for breaking if we don't throw anything new out)
        int thrown_out = 0; 
        
        // diagnostic. count the number of times OI is performed
        int count_oi = 0;

        // loop over observations
        for(int curr=0; curr < p; curr++) {

            if(debug) std::cout << "===> curr " << curr << " ===============" << std::endl;

            // -~- Can the observation be used as a centroid for the definition of inner and outer circles?

            // jump to next if observation not to be checked or already checked 
            if(obs_test[curr] != 1 || flags[curr] >= 0) {
                if(debug) std::cout << "..observation not suitable as a centroid " << curr << std::endl;
                continue;
            }

            // -~- Define outer/inner circles and which are the observations to test

            // get all neighbours that are close enough (inside outer circle)
            vec distances;
            ivec indices_global_outer_guess = tree.get_neighbours_with_distance(lats[curr], lons[curr], outer_radius, distances, true);
            if(indices_global_outer_guess.size() > num_max_outer) {

                int N = indices_global_outer_guess.size();
                std::vector<std::pair<float,int> > pairs(N);
                for(int i = 0; i < indices_global_outer_guess.size(); i++) {
                    pairs[i] = std::pair<float, int>(distances[i], indices_global_outer_guess[i]);
                }
                std::sort(pairs.begin(), pairs.end(), titanlib::sort_pair_first<float,int>());
                distances.clear();
                indices_global_outer_guess.clear();
                distances.resize(num_max_outer);
                indices_global_outer_guess.resize(num_max_outer);
                for(int i = 0; i < num_max_outer; i++) {
                    distances[i] = pairs[i].first;
                    indices_global_outer_guess[i] = pairs[i].second;
                }
            }
            if(debug) {
                int p_dist = distances.size();
                std::cout << "p_dist  " << p_dist << std::endl;
            }

            // set all the indices linking the different levels

            ivec indices_global_outer, indices_global_test, indices_outer_inner, indices_outer_test, indices_inner_test;
            bool res = titanlib::set_indices( indices_global_outer_guess, obs_test, flags, distances, inner_radius, -1, indices_global_outer, indices_global_test, indices_outer_inner, indices_outer_test, indices_inner_test);

            int p_outer = indices_global_outer.size();
            int p_inner = indices_outer_inner.size();
            int p_test = indices_global_test.size();
            if(debug) std::cout << "p_outer inner test " << p_outer << " " << p_inner << " " << "" << p_test << std::endl;

            // -~- Check if there are enough observations in the outer/inner circles

            if(p_outer < num_min_outer) {
                flags[curr] = flag_not_tested_outer;
                if(debug) {
                    std::cout << "@@isolated (outer) " << curr << std::endl;
                }
                continue; 
            }

//            if( p_inner < 2) {
//                flags[curr] = flag_not_tested_inner;
//                if(debug) std::cout << "@@isolated (inner) " << curr << std::endl;
//                continue;
//            } 
            
            // -~- Vectors on the outer circle

            vec lons_outer   = titanlib::subset( lons, indices_global_outer);
            vec elevs_outer  = titanlib::subset( elevs, indices_global_outer);
            vec lats_outer   = titanlib::subset( lats, indices_global_outer);
            vec eps2_outer   = titanlib::subset( eps2, indices_global_outer);
            vec w_outer      = titanlib::subset( w, indices_global_outer);
            vec t_outer      = titanlib::subset( t, indices_global_outer);
            
            // Debug for indices
            if(debug) {
                ivec flags_outer;
                flags_outer.resize( p_outer, na);
                for(int i=0; i < p_outer; i++) {
                   int g = indices_global_outer[i];
                   flags_outer[i] = flags[g];
                }
                std::cout << "curr lats lons elevs obs:" << std::endl;
                std::cout << curr << " " << lats[curr] << " " << lons[curr] << " " << elevs[curr] << " " << values[curr] << std::endl;
                std::cout << "indices_global_outer - i g lats lons elevs obs w flags:" << std::endl;
                for(int i=0; i<p_outer; i++) {
                    int g = indices_global_outer[i];
                    std::cout << std::setprecision(6) << i << " " << g << " " << lats[g] << " " << lons[g] << " " << elevs[g] << " " << values[g] << " " << w[g] << " " << flags[g] << std::endl;
                }
                std::cout << "outer - i lats lons elevs obs flags:" << std::endl;
                for(int i=0; i<p_outer; i++) {
                    std::cout << std::setprecision(6) << i << " " << lats_outer[i] << " " << lons_outer[i] << " " << elevs_outer[i] << " " << w_outer[i] << " " << flags_outer[i] << std::endl;
                }
                std::cout << "indices_outer_inner - l i lats lons elevs obs flags:" << std::endl;
                for(int l=0; l<p_inner; l++) {
                    int i = indices_outer_inner[l];
                    std::cout << std::setprecision(6) << l << " " << i << " " << lats_outer[i] << " " << lons_outer[i] << " " << elevs_outer[i] << " " << w_outer[i] << " " << flags_outer[i] << std::endl;
                }
                std::cout << "indices_outer_test - l m lats lons elevs obs flags:" << std::endl;
                for(int m=0; m<p_test; m++) {
                    int i = indices_outer_test[m];
                    std::cout << std::setprecision(6) << m << " " << i << " " << lats_outer[i] << " " << lons_outer[i] << " " << elevs_outer[i] << " " << w_outer[i] << " " << flags_outer[i] << std::endl;
                }
                std::cout << "indices_inner_test - l m lats lons elevs obs flags:" << std::endl;
                for(int m=0; m<p_test; m++) {
                    int l = indices_inner_test[m];
                    int i = indices_outer_inner[l];
                    std::cout << std::setprecision(6) << m << " " << l << " " << lats_outer[i] << " " << lons_outer[i] << " " << elevs_outer[i] << " " << w_outer[i] << " " << flags_outer[i] << std::endl;
                }
            }


            // -~- If all observations in the outer circle belongs to the same class then all are good

            int p_outer_w1 = 0;
            for(int i=0; i<p_outer; i++) 
                if ( w_outer[i] == 1) p_outer_w1++; 

            if ( p_outer_w1 == p_outer ||  p_outer_w1 == 0) {
                for(int m=0; m<p_test; m++) {
                    int g = indices_global_test[m];
                    flags[g] = 0;
                    if (debug) std::cout << " shortcut - outer circle has all w=0 or all w=1 - index " << g << std::endl;
                }
                continue;
            }
 
            // -~- SCT_dual on a selection of observations to test in the inner circle
            bool set_flag0 = iteration > 0;
            res = sct_dual_core( lats_outer, lons_outer, elevs_outer, w_outer, t_outer, min_horizontal_scale, max_horizontal_scale, kth_closest_obs_horizontal_scale, vertical_scale, eps2_outer, indices_global_test, indices_outer_inner, indices_outer_test, indices_inner_test, debug, na, set_flag0, thrown_out, flags); 
            // thrown_out is reset every new iteration

            // problems during the matrix inversion
            if ( !res) {
                flags[curr] = 100;
                if(debug) std::cout << " oi - flags=100 - index " << curr << std::endl; 
                continue;
            }

            count_oi++; // one more OI

        }  // end loop over observations

        double e_time0 = titanlib::clock();
        if (debug) std::cout << "SCT_dual loop - Removing " << thrown_out << " observations. Number of OI " << count_oi << " time=" << e_time0 - s_time0 << " secs" << std::endl;
        if(thrown_out == 0) {
            if ( iteration == 0) set_all_good = true;
            if(iteration + 1 < num_iterations)
                if (debug) std::cout << "Stopping early after " << iteration + 1<< " iterations" << std::endl;
            break;
        }
    
    } // end of SCT_dual iterations

    /* enter if the first iteration has been completed without flagging any bad observation
       The loop flags all the observations as good ones */
    if ( set_all_good) {
        int count_zero = 0;
        for(int curr=0; curr < p; curr++) {
          if ( flags[curr] == na && obs_test[curr] == 1) {
            flags[curr] = 0;
            count_zero++;
          }
        }
        if(debug) std::cout << "Remaining " << count_zero << " observations set to \"good\" " << std::endl;
    }

    /* check on observations missing QC flags */

    if(debug) std::cout << " +++++ check on observations missing QC flags ++++++++++++++++" << std::endl;
    double s_time0 = titanlib::clock();
        
    // reset this number each loop (this is for breaking if we don't throw anything new out)
    int thrown_out = 0; 
        
    // diagnostic. count the number of times OI is performed
    int count_oi = 0;

    // loop over observations
    for(int curr=0; curr < p; curr++) {

        if(debug) std::cout << "===> curr " << curr << " ===============" << std::endl;

        // -~- Can the observation be used as a centroid for the definition of inner and outer circles?

        if(obs_test[curr] != 1 || flags[curr] >= 0) {
            if(debug) std::cout << "..skip " << curr << std::endl;
            continue;
        }

        // -~- Define outer/inner circles and which are the observations to test

        // get all neighbours that are close enough (inside outer circle)
        vec distances;
        ivec indices_global_outer_guess = tree.get_neighbours_with_distance(lats[curr], lons[curr], outer_radius, distances, true);
        if(indices_global_outer_guess.size() > num_max_outer) {

            int N = indices_global_outer_guess.size();
            std::vector<std::pair<float,int> > pairs(N);
            for(int i = 0; i < indices_global_outer_guess.size(); i++) {
                pairs[i] = std::pair<float, int>(distances[i], indices_global_outer_guess[i]);
            }
            std::sort(pairs.begin(), pairs.end(), titanlib::sort_pair_first<float,int>());
            distances.clear();
            indices_global_outer_guess.clear();
            distances.resize(num_max_outer);
            indices_global_outer_guess.resize(num_max_outer);
            for(int i = 0; i < num_max_outer; i++) {
                distances[i] = pairs[i].first;
                indices_global_outer_guess[i] = pairs[i].second;
            }
        }

        // set all the indices linking the different levels
        // note that the only observation to test is curr
        ivec indices_global_outer, indices_global_test, indices_outer_inner, indices_outer_test, indices_inner_test;
        bool res = titanlib::set_indices( indices_global_outer_guess, obs_test, flags, distances, inner_radius, curr, indices_global_outer, indices_global_test, indices_outer_inner, indices_outer_test, indices_inner_test);
            
        int p_outer = indices_global_outer.size();
        int p_inner = indices_outer_inner.size();
        int p_test = indices_global_test.size();
        if(debug) std::cout << "p_outer inner test " << p_outer << " " << p_inner << " " << p_test << std::endl;
            
        // -~- Check if there are enough observations in the outer/inner circles

        if(p_outer < num_min_outer) {
            flags[curr] = flag_not_tested_outer;
            if(debug){
                std::cout << "@@isolated (outer) " << curr << std::endl;
            }
            continue; 
        }

//        if( p_inner < 2) {
//            flags[curr] = flag_not_tested_inner;
//            if(debug) std::cout << "@@isolated (inner) " << curr << std::endl;
//            continue;
//        }
            
        // -~- Vectors on the outer circle

        vec lons_outer   = titanlib::subset( lons, indices_global_outer);
        vec elevs_outer  = titanlib::subset( elevs, indices_global_outer);
        vec lats_outer   = titanlib::subset( lats, indices_global_outer);
        vec eps2_outer   = titanlib::subset( eps2, indices_global_outer);
        vec w_outer      = titanlib::subset( w, indices_global_outer);
        vec t_outer      = titanlib::subset( t, indices_global_outer);

        // -~- If all observations in the outer circle belongs to the same class then all are good

        int p_outer_w1 = 0;
        for(int i=0; i<p_outer; i++) 
            if ( w_outer[i] == 1) p_outer_w1++; 

        if ( p_outer_w1 == p_outer ||  p_outer_w1 == 0) {
            int g = indices_global_test[0];
            flags[g] = 0;
            if (debug) std::cout << " shortcut - outer circle has all w=0 or all w=1 - index " << g << std::endl;
            continue;
        }

        // -~- SCT_dual on a selection of observations to test in the inner circle

        res = sct_dual_core( lats_outer, lons_outer, elevs_outer, w_outer, t_outer, min_horizontal_scale, max_horizontal_scale, kth_closest_obs_horizontal_scale, vertical_scale, eps2_outer, indices_global_test, indices_outer_inner, indices_outer_test, indices_inner_test, debug, na, true, thrown_out, flags); 

        // problems during the matrix inversion
        if ( !res) {
            flags[curr] = 100;
            if(debug) std::cout << " oi - flags=100 - index " << curr << std::endl; 
            continue;
        }

        count_oi++; // one more OI

    }  // end loop over observations

    double e_time0 = titanlib::clock();
    if (debug) std::cout << "QC missing - Removing " << thrown_out << " observations. Number of OI " << count_oi << " time=" << e_time0 - s_time0 << " secs" << std::endl;


    /* final check on the bad observations
       it may happen that a good observation is flagged as a bad one because of the order
       the SCT_dual has been done and the uncertainty made in the estimation of the background 
       (remember that bad observations may have been used to get the background).
       This final check is made only on bad observations and uses only good observations. */

    if(debug) std::cout << " +++++ final check on the bad observations ++++++++++++++++" << std::endl;
    s_time0 = titanlib::clock();
        
    // reset this number each loop (this is for breaking if we don't throw anything new out)
    thrown_out = 0; 
        
    // diagnostic. count the number of times OI is performed
    count_oi = 0;

    // loop over observations
    for(int curr=0; curr < p; curr++) {

        if(debug) std::cout << "===> curr " << curr << " ===============" << std::endl;

        // -~- Can the observation be used as a centroid for the definition of inner and outer circles?
        if(obs_test[curr] != 1 || flags[curr] != 1) {
            if(debug) std::cout << "..skip " << curr << std::endl;
            continue;
        }

        // -~- Define outer/inner circles and which are the observations to test

        // get all neighbours that are close enough (inside outer circle)
        vec distances;
        ivec indices_global_outer_guess = tree.get_neighbours_with_distance(lats[curr], lons[curr], outer_radius, distances, true);
        if(indices_global_outer_guess.size() > num_max_outer) {

            int N = indices_global_outer_guess.size();
            std::vector<std::pair<float,int> > pairs(N);
            for(int i = 0; i < indices_global_outer_guess.size(); i++) {
                pairs[i] = std::pair<float, int>(distances[i], indices_global_outer_guess[i]);
            }
            std::sort(pairs.begin(), pairs.end(), titanlib::sort_pair_first<float,int>());
            distances.clear();
            indices_global_outer_guess.clear();
            distances.resize(num_max_outer);
            indices_global_outer_guess.resize(num_max_outer);
            for(int i = 0; i < num_max_outer; i++) {
                distances[i] = pairs[i].first;
                indices_global_outer_guess[i] = pairs[i].second;
            }
        }

        // set all the indices linking the different levels

        ivec indices_global_outer, indices_global_test, indices_outer_inner, indices_outer_test, indices_inner_test;
        bool res = titanlib::set_indices( indices_global_outer_guess, obs_test, flags, distances, inner_radius, curr, indices_global_outer, indices_global_test, indices_outer_inner, indices_outer_test, indices_inner_test);
            
        int p_outer = indices_global_outer.size();
        int p_inner = indices_outer_inner.size();
        int p_test = indices_global_test.size();
        if(debug) std::cout << "p_outer inner test " << p_outer << " " << p_inner << " " << p_test << std::endl;
            
        // -~- Decide if there are enough observations in the outer/inner circles

        if(p_outer < num_min_outer) {
            flags[curr] = flag_not_tested_outer;
            if(debug){
                std::cout << "@@isolated (outer) " << curr << std::endl;
            }
            continue; 
        }

//        if( p_inner < 2) {
//            flags[curr] = 11;
//            if(debug) std::cout << "@@isolated (inner) " << curr << std::endl;
//            continue;
//        }
            
        // -~- Set vectors on the outer circle

        vec lons_outer   = titanlib::subset( lons, indices_global_outer);
        vec elevs_outer  = titanlib::subset( elevs, indices_global_outer);
        vec lats_outer   = titanlib::subset( lats, indices_global_outer);
        vec eps2_outer   = titanlib::subset( eps2, indices_global_outer);
        vec w_outer      = titanlib::subset( w, indices_global_outer);
        vec t_outer      = titanlib::subset( t, indices_global_outer);

        // Debug for indices
        if(debug) {
            ivec flags_outer;
            flags_outer.resize( p_outer, na);
            for(int i=0; i < p_outer; i++) {
               int g = indices_global_outer[i];
               flags_outer[i] = flags[g];
            }
            std::cout << "curr lats lons elevs obs:" << std::endl;
            std::cout << curr << " " << lats[curr] << " " << lons[curr] << " " << elevs[curr] << " " << values[curr] << std::endl;
            std::cout << "indices_global_outer - i g lats lons elevs obs w flags:" << std::endl;
            for(int i=0; i<p_outer; i++) {
                int g = indices_global_outer[i];
                std::cout << std::setprecision(6) << i << " " << g << " " << lats[g] << " " << lons[g] << " " << elevs[g] << " " << values[g] << " " << w[g] << " " << flags[g] << std::endl;
            }
            std::cout << "outer - i lats lons elevs obs flags:" << std::endl;
            for(int i=0; i<p_outer; i++) {
                std::cout << std::setprecision(6) << i << " " << lats_outer[i] << " " << lons_outer[i] << " " << elevs_outer[i] << " " << w_outer[i] << " " << flags_outer[i] << std::endl;
            }
            std::cout << "indices_outer_inner - l i lats lons elevs obs flags:" << std::endl;
            for(int l=0; l<p_inner; l++) {
                int i = indices_outer_inner[l];
                std::cout << std::setprecision(6) << l << " " << i << " " << lats_outer[i] << " " << lons_outer[i] << " " << elevs_outer[i] << " " << w_outer[i] << " " << flags_outer[i] << std::endl;
            }
            std::cout << "indices_outer_test - l m lats lons elevs obs flag:" << std::endl;
            for(int m=0; m<p_test; m++) {
                int i = indices_outer_test[m];
                std::cout << std::setprecision(6) << m << " " << i << " " << lats_outer[i] << " " << lons_outer[i] << " " << elevs_outer[i] << " " << w_outer[i] << " " << flags_outer[i] << std::endl;
            }
            std::cout << "indices_inner_test - l m lats lons elevs obs flags:" << std::endl;
            for(int m=0; m<p_test; m++) {
                int l = indices_inner_test[m];
                int i = indices_outer_inner[l];
                std::cout << std::setprecision(6) << m << " " << l << " " << lats_outer[i] << " " << lons_outer[i] << " " << elevs_outer[i] << " " << w_outer[i] << " " << flags_outer[i] << std::endl;
            }
        }

        // -~- If all observations in the outer circle belongs to the same class then all are good

        int p_outer_w1 = 0;
        for(int i=0; i<p_outer; i++) 
            if ( w_outer[i] == 1) p_outer_w1++; 

        if ( p_outer_w1 == p_outer ||  p_outer_w1 == 0) {
            int g = indices_global_test[0];
            flags[g] = 0;
            if (debug) std::cout << " shortcut - outer circle has all w=0 or all w=1 - index " << g << std::endl;
            continue;
        }

        // -~- SCT_dual within the inner circle
        // note that the only observation to test is curr
        res = sct_dual_core( lats_outer, lons_outer, elevs_outer, w_outer, t_outer, min_horizontal_scale, max_horizontal_scale, kth_closest_obs_horizontal_scale, vertical_scale, eps2_outer, indices_global_test, indices_outer_inner, indices_outer_test, indices_inner_test, debug, na, true, thrown_out, flags); 

        // problems during the matrix inversion
        if ( !res) {
            flags[curr] = 100;
            if(debug) std::cout << " oi - flags=100 - index " << curr << std::endl; 
            continue;
        }

        count_oi++; // one more OI

    }  // end loop over observations

    e_time0 = titanlib::clock();
    if (debug) std::cout << "Re-check bad obs - Removing " << thrown_out << " observations. Number of OI " << count_oi << "time = " << e_time0 - s_time0 << " secs" << std::endl;

    //
    if(debug) {
        int count_bad = 0;
        int count_good = 0;
        int count_missing = 0;
        int count_iso_inner = 0;
        int count_iso_outer = 0;
        int count_fail_matinv = 0;
        int count_impossible = 0;
        // loop over observations
        for(int curr=0; curr < p; curr++) {
            if( flags[curr] == 1) {
                count_bad++;
            } else if( flags[curr] == 0) {
                count_good++;
            } else if ( flags[curr] == na) {
                count_missing++;
            } else if( accept_not_tested && (flags[curr] == flag_not_tested_inner)) {
                count_iso_inner++;
            } else if( accept_not_tested && (flags[curr] == flag_not_tested_outer))  {
                count_iso_outer++;
            } else if( flags[curr] == 100) {
                count_fail_matinv++;
            } else {
                count_impossible++;
            }
        }
        if(accept_not_tested)
            std::cout << std::setprecision(3) << "summary - # TOT good bad missing: " << p << " " << count_good << " " << count_bad << " " << count_missing << std::endl; 
        else
            std::cout << std::setprecision(3) << "summary - # TOT good bad missing isolated(inner) isolated(outer): " << p << " " << count_good << " " << count_bad << " " << count_missing << " " << count_iso_inner << " " << count_iso_outer << std::endl;
        if(count_fail_matinv > 0) 
            std::cout << std::setprecision(3) << "!!!! failure in matrix inversion: " << count_fail_matinv << std::endl; 
        if(count_impossible > 0) 
            std::cout << std::setprecision(3) << "!!!! unknown flag: " << count_impossible << std::endl; 
    }
    //
    
    if(debug) std::cout << ">> Total Time " << e_time0 - s_time << "secs" << std::endl;

    //
    return flags;
}
// end SCT_dual //


//----------------------------------------------------------------------------//
// SCT dual WITHIN THE INNER CIRCLE

// Smallets SCT_dual unit
bool sct_dual_core( const vec& lats, 
                    const vec& lons, 
                    const vec& elevs, 
                    const vec& w, 
                    const vec& t, 
                    float Dh_min, 
                    float Dh_max, 
                    int kth_close, 
                    float Dz, 
                    const vec& eps2, 
                    const ivec& indices_global_test, 
                    const ivec& indices_outer_inner, 
                    const ivec& indices_outer_test, 
                    const ivec& indices_inner_test, 
                    bool debug, 
                    float na, 
                    bool set_flag0, 
                    int& thrown_out,
                    ivec& flags) {
/*-----------------------------------------------------------------------------
 Spatial consistency test for dichotomous (yes/no) variables to test within the
 inner circle, considering all the observations in the outer circle.

 The inner and outer circles are defined with respect to a centroid observation.

 Each observation to test is divided into "yes, the event happened" or 
 "no, the event did not happen" (w vector: yes=1, no=0)
 In the outer circle, there are two sub-networks of observations. 
 The yes-network (w1) and the no-network (w0).
 We assume that information is distributed in space through correlation functions, 
 analougously to what happen for statistical interpolation with error correlation
 matrices.
 At each observation location, the leave-one-out integral data influence (idiv)
 is obtained for each of the two sub-networks. The two values w1_idiv and w0_idiv 
 are then obtained. The larger of those two values corresponds to the best guess
 between "yes" and "no" derived by the neighbours, which is also a guess that is 
 independent from the actual observation. If there is a mismatch between the 
 actual obtained w value and the best guess, then we have a canditate bad 
 observation.
 Note that when computing a best guess we assume all buddies are good ones. Then, 
 a bad buddy may cause a good observation to be misjudged as bad. For this reason
 we allow for the rejection of just one bad observation at a time. We reject the 
 one where the mismatch between actual value and best guess is the worst.

 To quantify the mismatches, we use the relative information content of 
 w1_idiv with respect to w0_idiv (and viceversa I(0wrt1)):
  I(1wrt0) = w1_idiv * log( w1_idiv / w0_idiv)
 (see e.g. page 12 of Tarantola's Book "Inverse problem theory" (2005).)
 which is equivalent to defining the distance between w1_idiv and w0_idiv not
 through difference but in a more complicated way. 
 Note than I(...) is called "z" in the code.
 If no mismatches are found, then all observations are good ones.
 
 In areas of transition between "yes" and "no" regions, a good observation can be 
 mistaken for a bad one. For this reason, the user is allowed to set a threshold "t"
 for the relative information content and an observation is flagged as bad only if
 the relative information content I is larger than the threshold.

 Note that if both w1_idiv and w0_idiv are small numbers (close to zero), then
 we are dealing with an observation that is rather distant from its buddies.
 The distance is specified with respect to the metric defined by the correlation
 functions.
 In this case, better to leave the judgment suspended and consider the 
 observation good.
 A small measure of optimism in ignorance is fashionable these days, 
 we can't help but adapt.

  Output values:

   flags, global vector. return updated quality control flag values (1 = bad, 0 = good)
   
   thrown_out, integer. update the number of observations flagged as bad ones 

-----------------------------------------------------------------------------*/

    using namespace boost::numeric::ublas;

    // init
    int p_outer = w.size();
    int p_test  = indices_global_test.size();

    /* observations where leave-one-out IDIs are both small are difficult to judge
       because they do not have close-enough buddies
       we skip them. This is the threshold */
    float w_idiv_min = 0.45;
 
    /* two more indices 
      ... w0[i] = na when w[i] = 1;  w0[i] = j for the j-th element of w where w = 0
      ... w1[i] = na when w[i] = 0;  w1[i] = k for the k-th element of w where w = 1
    */
    ivec indices_outer_w0( p_outer, na);
    ivec indices_outer_w1( p_outer, na);
    int j0=0;
    int j1=0;
    for(int i=0; i < p_outer; i++) {
      if ( w[i] == 0) {
          indices_outer_w0[i] = j0;
          j0++;
      }
      if ( w[i] == 1) {
          indices_outer_w1[i] = j1;
          j1++;
      }
    }
    int p_outer_w0 = j0;
    int p_outer_w1 = j1;
    if (debug) {
        std::cout << " p_outer_w0 p_outer_w1 = " << p_outer_w0 << " " << p_outer_w1 << std::endl;
        for(int i=0; i < p_outer; i++) {
            std::cout << " i indices_outer_w0 indices_outer_w1 " << indices_outer_w0[i] << " " << indices_outer_w1[i] << std::endl;
        }
    }

    /* Compute Dh. 
     The location-dependent horizontal de-correlation lenght scale 
     used for the background error correlation matrix */
    boost::numeric::ublas::matrix<float> disth(p_outer, p_outer);
    boost::numeric::ublas::matrix<float> distz(p_outer, p_outer);
    boost::numeric::ublas::vector<float> Dh(p_outer);
//    boost::numeric::ublas::matrix<float> disth_w0(p_outer_w0, p_outer_w0);
//    boost::numeric::ublas::matrix<float> distz_w0(p_outer_w0, p_outer_w0);
//    boost::numeric::ublas::matrix<float> disth_w1(p_outer_w1, p_outer_w1);
//    boost::numeric::ublas::matrix<float> distz_w1(p_outer_w1, p_outer_w1);

    for(int i=0; i < p_outer; i++) {
        vec Dh_vector(p_outer);
//        int i0 = indices_outer_w0[i];
//        int i1 = indices_outer_w1[i];
        int k = 0;
        for(int j=i; j < p_outer; j++) {
//            int j0 = indices_outer_w0[j];
//            int j1 = indices_outer_w1[j];
            disth(i, j) = titanlib::calc_distance( lats[i], lons[i], lats[j], lons[j]);
            distz(i, j) = fabs( elevs[i] - elevs[j]);
//            if ( i1 != na && j1 != na) {
//                disth_w1(i1, j1) = disth(i, j);
//                distz_w1(i1, j1) = distz(i, j);
//            } else if ( i0 != na && j0 != na) {
//                disth_w0(i0, j0) = disth(i, j);
//                distz_w0(i0, j0) = distz(i, j);
//            }
            if(i != j) {
                disth(j, i) = disth(i, j);
                distz(j, i) = distz(i, j);
//                if ( i1 != na && j1 != na) {
//                    disth_w1(j1, i1) = disth(i, j);
//                    distz_w1(j1, i1) = distz(i, j);
//                } else if ( i0 != na && j0 != na) {
//                    disth_w0(j0, i0) = disth(i, j);
//                    distz_w0(j0, i0) = distz(i, j);
//                }
                Dh_vector[k] = disth(i, j);
                k++;
            }
        }
        for(int j=0; j <= i; j++) {
            Dh_vector[k] = disth(i, j);
            k++;
        }
        // find distance to the k-th closest observations
        Dh(i) = titanlib::find_k_closest(Dh_vector, kth_close);
        if(debug) std::cout << "i Dh " << i << " " << Dh(i) << std::endl;
    }

    float Dh_mean = std::accumulate(std::begin(Dh), std::end(Dh), 0.0) / Dh.size();
    if(Dh_mean < Dh_min) {
        if(debug) std::cout << "Dh_mean (<min_horizontal_scale) " << Dh_mean << std::endl;
        Dh_mean = Dh_min;
    }
    if(Dh_mean > Dh_max) {
        if(debug) std::cout << "Dh_mean (>max_horizontal_scale) " << Dh_mean << std::endl;
        Dh_mean = Dh_max;
    }
    if(debug) std::cout << "Dh_mean " << Dh_mean << std::endl;
    if(debug) std::cout << "Dz " << Dz << std::endl;
    
    /* Compute S + eps2*I and store it in S 
       S_w1 made by only points where w = 1
       S_w0 made by only points where w = 0
    */
    boost::numeric::ublas::matrix<float> S(p_outer,p_outer);
    boost::numeric::ublas::matrix<float> S_w0(p_outer_w0,p_outer_w0);
    boost::numeric::ublas::matrix<float> S_w1(p_outer_w1,p_outer_w1);
    boost::numeric::ublas::matrix<float> Sinv_w0(p_outer_w0,p_outer_w0);
    boost::numeric::ublas::matrix<float> Sinv_w1(p_outer_w1,p_outer_w1);
    for(int i=0; i < p_outer; i++) {
        int i0 = indices_outer_w0[i];
        int i1 = indices_outer_w1[i];
        for(int j=i; j < p_outer; j++) {
            int j0 = indices_outer_w0[j];
            int j1 = indices_outer_w1[j];
            float value = std::exp(-.5 * std::pow((disth(i, j) / Dh_mean), 2) - .5 * std::pow((distz(i, j) / Dz), 2));
            if(i==j) { // weight the diagonal, this also ensure an invertible matrix
                value = value + eps2[i];
            } else {
                S(j, i) = value;
            }
            S(i, j) = value;
            if ( i1 != na && j1 != na) {
                S_w1(i1, j1) = value;
                if ( i1 != j1) S_w1(j1, i1) = value;
            } else if ( i0 != na && j0 != na) {
                S_w0(i0, j0) = value;
                if ( i0 != j0) S_w0(j0, i0) = value;
            }
        }
    }
    
    // Invert S matrices
    bool b = invert_matrix( S_w0, Sinv_w0);
    if( !b)  return( false);
    b = invert_matrix( S_w1, Sinv_w1);
    if( !b)  return( false);

    // Definitions

    boost::numeric::ublas::vector<float> Sinv_d_w0(p_outer_w0);
    boost::numeric::ublas::vector<float> Sinv_d_w1(p_outer_w1);

    // re-obtains S and sum-up the rows of the Inv(S_w) matrices

    for(int i=0; i<p_outer; i++) {
        S(i,i) -= eps2[i];
        int i0 = indices_outer_w0[i];
        int i1 = indices_outer_w1[i];
        if ( i1 != na) {
            float acc = 0;
            for(int j1=0; j1<p_outer_w1; j1++)
                acc += Sinv_w1(i1,j1);
            Sinv_d_w1(i1) = acc;
        } else {
            float acc = 0;
            for(int j0=0; j0<p_outer_w0; j0++)
                acc += Sinv_w0(i0,j0);
            Sinv_d_w0(i0) = acc;
        }
    }

    // Run the real check over the observations to test

    float zmx = na;
    int mmx = na;
    if(debug) std::cout << "z_test - i lon lat w w1_idiv w0_idiv" << std::endl;
    for(int m=0; m<p_test; m++) {
        int i = indices_outer_test[m];
        int i0 = indices_outer_w0[i];
        int i1 = indices_outer_w1[i];
        // compute leave-one-out IDIs, one considering w=0 only and the other with w=1 only
        float w0_idiv = 0; 
        float w1_idiv = 0; 
        for(int j=0; j<p_outer; j++){
            int j0 = indices_outer_w0[j];
            int j1 = indices_outer_w1[j];
            if ( i1 != na && j0 != na) {
                w0_idiv += S(i, j) * Sinv_d_w0(j0);
            } else if ( i0 != na && j1 != na) {
                w1_idiv += S(i, j) * Sinv_d_w1(j1); 
            }
        }
        if ( i1 != na) {
          w1_idiv = 1. - Sinv_d_w1(i1) / Sinv_w1(i1, i1);
//          if (w1_idiv < w_idiv_min) w1_idiv = w_idiv_min;
          if (w1_idiv <= 0) w1_idiv = 0.001;
        } else {
          w0_idiv = 1. - Sinv_d_w0(i0) / Sinv_w0(i0, i0);
//          if (w0_idiv < w_idiv_min) w0_idiv = w_idiv_min;
          if (w0_idiv <= 0) w0_idiv = 0.001;
        }
        float z = na;
        float z0wrt1 = w0_idiv * log( w0_idiv / w1_idiv);
        float z1wrt0 = w1_idiv * log( w1_idiv / w0_idiv);
        // verify that the observation has buddies nearby
        if ( w1_idiv >= w_idiv_min || w0_idiv >= w_idiv_min) {
            // candidate bad if w=1 but the neighbours consider more likely w=0
            if ( w[i] == 1 && w0_idiv > w1_idiv && z0wrt1 > t[i]) {
                z = z0wrt1;
            // candidate bad if w=0 but the neighbours consider more likely w=1
            } else if ( w[i] == 0 && w1_idiv > w0_idiv && z1wrt0 > t[i]) {
                z = z1wrt0;
            }
            // real bad is the observation with larger z
            if ( z != na)  {
                if ( zmx == na) {
                    zmx = z;
                    mmx = m;
                } else if ( z > zmx ) {
                    zmx = z;
                    mmx = m;
                }
            }
        }
        if(debug) std::cout << i << " " << w[i] << " " << lons[i] << " " << lats[i] << " " << w1_idiv << " " << w0_idiv << " " << z1wrt0 << " " << z0wrt1 << " " << z << " " << zmx << " " << mmx << std::endl;
    }
    
    // Take a decision

    if ( mmx != na) {
        // flag just the worse z 
        int g = indices_global_test[mmx];
        int i = indices_outer_test[mmx];
        flags[g] = 1;
        thrown_out++;
        if(debug) std::cout << std::setprecision(6) << "SCT_dual failed - flag as bad - index w z " << g << " " << w[i] << " " << zmx << std::endl;
    // ... BUT if no mismatch between obtained and expected w's occur, then flag ALL as good
    } else if (set_flag0) {
        for(int m=0; m<p_test; m++) {
            int g = indices_global_test[m];
            flags[g] = 0;
            if(debug) std::cout << "all set to good - g " << g << std::endl;
        }
    }
    
    // normal exit
    return true;       

}
