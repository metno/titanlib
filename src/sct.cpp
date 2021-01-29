#include <vector>
#include <math.h>
#include "titanlib.h"
#include <assert.h>
#include <iostream>
#include <numeric>
#include <exception>
#include <gsl/gsl_multimin.h>
#include <gsl/gsl_linalg.h>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/vector_proxy.hpp>
#include <boost/numeric/ublas/triangular.hpp>
#include <boost/numeric/ublas/lu.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/ublas/matrix.hpp>


//=============================================================================

// helpers
bool invert_matrix (const boost::numeric::ublas::matrix<float>& input, boost::numeric::ublas::matrix<float>& inverse);

bool set_indices( const ivec& indices_global_outer_guess, const ivec& obs_test, const ivec& flags, const vec& dist_outer_guess, float inner_radius, int test_just_this, ivec& indices_global_outer, ivec& indices_global_test, ivec& indices_outer_inner, ivec& indices_outer_test, ivec& indices_inner_test);

// SCT within the inner circle

bool sct_core( const vec& lats, const vec& lons, const vec& elevs, const vec& yo, const vec& yb, float Dh_min, float Dh_max, int kth_close, float Dz, const vec& eps2, float minp, float maxp, const vec& mina, const vec& maxa, const vec& minv, const vec& maxv, const vec& tpos, const vec& tneg, const ivec& indices_global_test, const ivec& indices_outer_inner, const ivec& indices_outer_test, const ivec& indices_inner_test, bool debug, float na, bool set_flag0, int& thrown_out, vec& scores, ivec& flags);

// background (first-guess estimates) elaboration

vec background( std::string background_elab_type, const vec& elevs, const vec& values, int num_min_prof, float min_elev_diff, float value_minp, float value_maxp, const vec& external_background_values, const ivec& indices_global_outer, bool debug);

vec compute_vertical_profile(const vec& elevs, const vec& oelevs, const vec& values, int num_min_prof, double min_elev_diff, bool debug);

vec compute_vertical_profile_Theil_Sen(const vec& elevs, const vec& oelevs, const vec& values, int num_min_prof, double min_elev_diff, bool debug);

double basic_vertical_profile_optimizer_function(const gsl_vector *v, void *data);

vec basic_vertical_profile(const int n, const double *elevs, const double t0, const double gamma);

double vertical_profile_optimizer_function(const gsl_vector *v, void *data);

vec vertical_profile(const int n, const double *elevs, const double t0, const double gamma, const double a, const double h0, const double h1i);

//=============================================================================

//+ SCT - Spatial Consistency Test
ivec titanlib::sct( const vec& lats,
                    const vec& lons,
                    const vec& elevs,
                    const vec& values,
                    const ivec& obs_to_check,
                    const vec& background_values,
                    std::string background_elab_type,
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
                    const vec& values_mina,
                    const vec& values_maxa,
                    const vec& values_minv,
                    const vec& values_maxv,
                    const vec& eps2,
                    const vec& tpos,
                    const vec& tneg,
                    bool debug,
                    vec& scores) {
/*

 -+-+ Spatial Consistency Test +-+-

 + Description: 
    Flag observations that are (likely) affected by gross measurement errors (GE) based on neighbouring observations. Observations are affected by GEs if their values are (a) [not related to the actual atmospheric state] OR (b) [affected by such large representativeness errors (REs) that they are difficult to reconstruct using the neighbouring observations].
  For simplicity, we will refer to observations affected by GE as bad observations, viceversa good observations are those not affected by GEs. 

 + Reference: 
    Lussana, C., Uboldi, F. and Salvati, M.R. (2010), A spatial consistency test for surface observations from mesoscale meteorological networks. Q.J.R. Meteorol. Soc., 136: 1075-1088. doi:10.1002/qj.622 (Reference Abbreviation: LUS10)
 See also the wiki-pages https://github.com/metno/titanlib/wiki/Spatial-consistency-test

 + Algorithm: 
    The big task of applying the SCT over the whole observational network is splitted into several smaller tasks.
  The SCT is applied many times in sequence. 
  Each time, the SCT is applied over a small region centered on an observation location (i.e. the centroid observation).
  Based on the centroid, we define an inner circle and an outer circle.
  The outer circle is used to obtain a spatial trend of the observed variable over the region, the so-called background.
  There are many ways to compute such background values.
  The observations to test are those within the inner circle that have not been tested yet.
  The quality of those observations is determined by comparing them against estimates calculated considering all observations within the outer circle. Considering an observations, we use an estimate that depends on the observed value (analysis) and one that is indipendent of it (leave-one-out cross-validation analysis, or cv-analysis).
  See also the comments of the function "sct_core" for more details.
  
  The application of the sequence of SCTs over the whole observational network is repeated several times, until no GEs are found within the observations tested. 
  In fact, the core of SCT is to test a single observation against an independent estimate that is obtained assuming all the neighbouring observations are good ones. As a consequence, a bad observation can influence the results of some tests on neighbouring observations before being rejected.
  We apply special care to avoid "false positives" (good flagged as bad).
  For each SCT over a centroid, we flag as bad only the observation showing the worst result.
  In addition, an extra-loop is done over the bad observations where only the good ones are considered. We may "save" bad observations.
  Analogously, to avoid "false negatives" (bad flagged as good), the last loop is done over the bad observations where only the good ones are considered. We may "save" bad observations during this last loop.

  (*) Notes on prior information (see also comments in sct_core): the user is expected to provide two different range of values, that inform about the desired tolerance level of the SCT.
    i.  The range of admissible values (values_mina,values_maxa). A range for each observation. The observed value have to be within the corresponding range (e.g. observation=2, then OK if values_mina=0 & values_maxa=4; not OK if values_mina=10 & values_maxa=14). A cv-analysis outside this range warns us about possible problems with the observation.
        a.   The range of plausible values (value_minp,value_maxp) is derived from the range of admissible values and it is used to check the analyses and backgrounds. Outside this range the values do not make sense.
    ii. The range of valid values (values_minv,values_maxv). A range for each observation. The observed value have to be within the corresponding range. If the cv-analysis is within this range, then we assume the corresponding observation is a good one, regardless of the SCT score.
  Examples:
    Temperature: 3 observed values(degC)=( 10, 15, 5). (values_mina,values_maxa) = observed values plus/minus 15 degC. (values_minv,values_maxv) = observed values plus/minus 1 degC.
    Precipitation: 3 observed values(mm)=( 0, 50, 1). (values_mina,values_maxa) = observed values plus/minus either 10 mm or 50% of the value. (values_minv,values_maxv) = observed values plus/minus either 1 mm or 10% of the value.

  (*) Notes on the data structures (see also comments in set_indices): we have defined a hierarchy of data structures
    global -> outer -> inner -> test
   global vectors are the input vectors. The other data structures refer to a specific iteration when a centroid observation has been identified. outer are the vectors of data on the outer circle, not yet flagged as bad. inner are the data on the inner circle, not yet flagged as bad. test data are the observations in the inner circle and still to test.
   We use the same counter variables for the data structures: g -> global; i,j -> outer; l -> inner; m -> test
   There are vectors of indices matching the different structures. 
   For example, indices_outer_inner matches the outer and the inner vectors. indices_outer_inner is a vector of the dimension of the number of observations in the inner circle. The instruction i=indices_outer_inner[l] maps the l-th element of inner-type vectors as the i-th element of outer-type vectors.

 (*) Notes on the thresholds. Different thresholds can be used for different observations. The thresholds vectors have the same dimensions of the observations vectors (global). For each observation, two thresholds (tpos,tneg) have to be set by the user. tpos is the threshold used when the observed value is higher than the cv-analysis, otherwise tneg is used. 

 + Pseudo-code:
 
  Loop over the number of iterations prescribed
    Loop over the p observations (observation index is "curr"). 
      Decide if the observation is suitable for a centroid
      Define outer/inner cirles, observations to test and set all the vector of indices
      Check if the observation is isolated
      Compute the background
      If the deviations between observations and background is small enough that further steps are not needed and the observations to test are flagged as good ones
      SCT over the inner circle, either flag the worse observation as bad or flag all of them as good ones
      NOTE: observations can be flagged as good ones only after the first iteration (after rejecting the worst)
  
  Loop over the observations without QC flags (nor good neither bad) and use each of them as a centroid

  Loop over the bad observations and use each of them as a centroid
  
 + Returned values:

  scores. only for observations rejected by the SCT, z-score (See sct_core)

  flags. -999 = not checked; 0 = passed (good); 1 = failed (bad); 11 = isolated (<2 inside inner); 12 = isolated (<num_min_outer inside outer)

*/
    double s_time = titanlib::util::clock();
    
    const int p = values.size();
    if( lats.size() != p || lons.size() != p || elevs.size() != p || values.size() != p || tpos.size() != p || tneg.size() != p || eps2.size() != p || values_mina.size() != p || values_maxa.size() != p || values_minv.size() != p || values_maxv.size() != p)
        throw std::runtime_error("Dimension mismatch");
    if(num_min_outer < 2)
        throw std::invalid_argument("num_min_outer must be > 1");
    if(num_max_outer < num_min_outer)
        throw std::invalid_argument("num_max_outer must be > num_min_outer");
    if(num_iterations < 1)
        throw std::invalid_argument("num_iterations must be >= 1");
    if(min_elev_diff <= 0)
        throw std::invalid_argument("min_elev_diff must be > 0");
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
    if( background_elab_type != "vertical_profile" && 
        background_elab_type != "vertical_profile_Theil_Sen" && 
        background_elab_type != "mean_outer_circle" && 
        background_elab_type != "median_outer_circle" &&
        background_elab_type != "external")
        throw std::invalid_argument("background_elab_type must be one of vertical_profile, vertical_profile_Theil_Sen, mean_outer_circle, median_outer_circle or external");
    if(background_elab_type == "vertical_profile" && num_min_prof<0)
        throw std::invalid_argument("num_min_prof must be >=0");
    if(background_elab_type == "external" &&  background_values.size() != p)
        throw std::runtime_error("Background vector dimension mismatch");

    // initializations
    float na = -999.; // code for "Not Available". Any type of missing data
    ivec flags( p, na);
    ivec obs_test( p, 1);

    bool set_all_good = false;
    
    scores.clear();
    scores.resize( p, na);
    
    if (obs_to_check.size() == p) {
        for(int g=0; g < p; g++) 
          obs_test[g] = obs_to_check[g];
    }

    /* plausible values are inferred from admissible values
       note: a range check should be performed before doing the SCT.
             The SCT can mistake obvious large errors for good observations 
             if they happen to be close enough that they support each other */
    float value_minp = values_mina[0];
    float value_maxp = values_maxa[0];
    for(int g=1; g < p; g++) {
        if ( values_mina[g] < value_minp) value_minp = values_mina[g];
        if ( values_maxa[g] > value_maxp) value_maxp = values_maxa[g];
    }
    if(debug) std::cout << "=============================================== " << p << std::endl;
    if(debug) std::cout << "Number of observations to test is " << p << std::endl;
    if(debug) std::cout << "Min and max acceptable values are " << value_minp << " " << value_maxp << std::endl;

    // KDtree has to do with fast computation of distances
    titanlib::KDTree tree(lats, lons);

    //
    if(debug) {
        std::cout << "g lats lons obs:" << std::endl;
        for(int g=0; g < p; g++) {
            std::cout << g << " " << lats[g] << " " << lons[g] << " " << values[g] << std::endl;
        }
    }

    // SCT iterations
    for(int iteration = 0; iteration < num_iterations; iteration++) {

        if(debug) std::cout << " +++++ Iteration " << iteration << " ++++++++++++++++" << std::endl;

        double s_time0 = titanlib::util::clock();
        
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
            ivec indices_global_outer_guess = tree.get_neighbours_with_distance(lats[curr], lons[curr], outer_radius, num_max_outer, true, distances);
            if(debug) {
                int p_dist = distances.size();
                std::cout << "p_dist  " << p_dist << std::endl;
            }

            // set all the indices linking the different levels

            ivec indices_global_outer, indices_global_test, indices_outer_inner, indices_outer_test, indices_inner_test;
            bool res = set_indices( indices_global_outer_guess, obs_test, flags, distances, inner_radius, -1, indices_global_outer, indices_global_test, indices_outer_inner, indices_outer_test, indices_inner_test);

            int p_outer = indices_global_outer.size();
            int p_inner = indices_outer_inner.size();
            int p_test = indices_global_test.size();
            if(debug) std::cout << "p_outer inner test " << p_outer << " " << p_inner << " " << "" << p_test << std::endl;

            // -~- Check if there are enough observations in the outer/inner circles

            if(p_outer < num_min_outer) {
                flags[curr] = 12;
                if(debug) std::cout << "@@isolated (outer) " << curr << std::endl;
                continue; 
            }

            if( p_inner < 2) {
                flags[curr] = 11;
                if(debug) std::cout << "@@isolated (inner) " << curr << std::endl;
                continue;
            }
            
            // -~- Vectors on the outer circle

            vec lons_outer   = titanlib::util::subset(lons, indices_global_outer);
            vec elevs_outer  = titanlib::util::subset(elevs, indices_global_outer);
            vec lats_outer   = titanlib::util::subset(lats, indices_global_outer);
            vec values_outer = titanlib::util::subset(values, indices_global_outer);
            vec eps2_outer   = titanlib::util::subset(eps2, indices_global_outer);
            vec tpos_outer   = titanlib::util::subset(tpos, indices_global_outer);
            vec tneg_outer   = titanlib::util::subset(tneg, indices_global_outer);
            vec mina_outer   = titanlib::util::subset(values_mina, indices_global_outer);
            vec minv_outer   = titanlib::util::subset(values_minv, indices_global_outer);
            vec maxa_outer   = titanlib::util::subset(values_maxa, indices_global_outer);
            vec maxv_outer   = titanlib::util::subset(values_maxv, indices_global_outer);
            
            // Debug for indices
            if(debug) {
                ivec flags_outer;
                flags_outer.resize( p_outer, na);
                for(int i=0; i < p_outer; i++) {
                   int g = indices_global_outer[i];
                   flags_outer[i] = flags[g];
                }
                std::cout << "curr lats lons obs:" << std::endl;
                std::cout << curr << " " << lats[curr] << " " << lons[curr] << " " << values[curr] << std::endl;
                std::cout << "indices_global_outer - i g lats lons obs:" << std::endl;
                for(int i=0; i<p_outer; i++) {
                    int g = indices_global_outer[i];
                    std::cout << std::setprecision(6) << i << " " << g << " " << lats[g] << " " << lons[g] << " " << values[g] << " " << flags[g] << std::endl;
                }
                std::cout << "outer - i lats lons obs:" << std::endl;
                for(int i=0; i<p_outer; i++) {
                    std::cout << std::setprecision(6) << i << " " << lats_outer[i] << " " << lons_outer[i] << " " << values_outer[i] << " " << flags_outer[i] << std::endl;
                }
                std::cout << "indices_outer_inner - l i lats lons obs:" << std::endl;
                for(int l=0; l<p_inner; l++) {
                    int i = indices_outer_inner[l];
                    std::cout << std::setprecision(6) << l << " " << i << " " << lats_outer[i] << " " << lons_outer[i] << " " << values_outer[i] << " " << flags_outer[i] << std::endl;
                }
                std::cout << "indices_outer_test - l m lats lons obs:" << std::endl;
                for(int m=0; m<p_test; m++) {
                    int i = indices_outer_test[m];
                    std::cout << std::setprecision(6) << m << " " << i << " " << lats_outer[i] << " " << lons_outer[i] << " " << values_outer[i] << " " << flags_outer[i] << std::endl;
                }
                std::cout << "indices_inner_test - l m lats lons obs:" << std::endl;
                for(int m=0; m<p_test; m++) {
                    int l = indices_inner_test[m];
                    int i = indices_outer_inner[l];
                    std::cout << std::setprecision(6) << m << " " << l << " " << lats_outer[i] << " " << lons_outer[i] << " " << values_outer[i] << " " << flags_outer[i] << std::endl;
                }
            }


            // -~- Compute the background 
            vec bvalues_outer = background( background_elab_type, elevs_outer, values_outer, num_min_prof, min_elev_diff, value_minp, value_maxp, background_values, indices_global_outer, debug);

            if(debug) std::cout << "... background ok ..." << std::endl;

            // -~- If deviations between backgrounds and observations are small then flag = 0

            /* if observations and background are almost identical (within the range of valid estimates),
               then take a shortcut and flag = 0 for all observations in the inner circle */

            bool small_innov = true;
            for(int m=0; m<p_test; m++) {
                int j = indices_outer_test[m];
                if(bvalues_outer[j] < minv_outer[j] || bvalues_outer[j] > maxv_outer[j]) 
                    small_innov = false;
            }

            if (small_innov) {
                for(int m=0; m<p_test; m++) {
                    int g = indices_global_test[m];
                    flags[g] = 0;
                    if (debug) std::cout << " small_innov - index " << g << std::endl;
                }
                continue;
            }
 
            // -~- SCT on a selection of observations to test in the inner circle
            bool set_flag0 = iteration > 0;
            res = sct_core( lats_outer, lons_outer, elevs_outer, values_outer, bvalues_outer, min_horizontal_scale, max_horizontal_scale, kth_closest_obs_horizontal_scale, vertical_scale, eps2_outer, value_minp, value_maxp, mina_outer, maxa_outer, minv_outer, maxv_outer, tpos_outer, tneg_outer, indices_global_test, indices_outer_inner, indices_outer_test, indices_inner_test, debug, na, set_flag0, thrown_out, scores, flags);
            // thrown_out is reset every new iteration

            // problems during the matrix inversion
            if ( !res) {
                flags[curr] = 100;
                if(debug) std::cout << " oi - flags=100 - index " << curr << std::endl; 
                continue;
            }

            count_oi++; // one more OI

        }  // end loop over observations

        std::cout << "SCT loop - Removing " << thrown_out << " observations. Number of OI " << count_oi << std::endl;
        double e_time0 = titanlib::util::clock();
        std::cout << e_time0 - s_time0 << " secs" << std::endl;
        if(thrown_out == 0) {
            if ( iteration == 0) set_all_good = true;
            if(iteration + 1 < num_iterations)
                std::cout << "Stopping early after " << iteration + 1<< " iterations" << std::endl;
            break;
        }
    
    } // end of SCT iterations

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
    double s_time0 = titanlib::util::clock();
        
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
        ivec indices_global_outer_guess = tree.get_neighbours_with_distance(lats[curr], lons[curr], outer_radius, num_max_outer, true, distances);

        // set all the indices linking the different levels
        // note that the only observation to test is curr
        ivec indices_global_outer, indices_global_test, indices_outer_inner, indices_outer_test, indices_inner_test;
        bool res = set_indices( indices_global_outer_guess, obs_test, flags, distances, inner_radius, curr, indices_global_outer, indices_global_test, indices_outer_inner, indices_outer_test, indices_inner_test);
            
        int p_outer = indices_global_outer.size();
        int p_inner = indices_outer_inner.size();
        int p_test = indices_global_test.size();
        if(debug) std::cout << "p_outer inner test " << p_outer << " " << p_inner << " " << p_test << std::endl;
            
        // -~- Check if there are enough observations in the outer/inner circles

        if(p_outer < num_min_outer) {
            flags[curr] = 12;
            if(debug) std::cout << "@@isolated (outer) " << curr << std::endl;
            continue; 
        }

        if( p_inner < 2) {
            flags[curr] = 11;
            if(debug) std::cout << "@@isolated (inner) " << curr << std::endl;
            continue;
        }
            
        // -~- Vectors on the outer circle

        vec lons_outer   = titanlib::util::subset(lons, indices_global_outer);
        vec elevs_outer  = titanlib::util::subset(elevs, indices_global_outer);
        vec lats_outer   = titanlib::util::subset(lats, indices_global_outer);
        vec values_outer = titanlib::util::subset(values, indices_global_outer);
        vec eps2_outer   = titanlib::util::subset(eps2, indices_global_outer);
        vec tpos_outer   = titanlib::util::subset(tpos, indices_global_outer);
        vec tneg_outer   = titanlib::util::subset(tneg, indices_global_outer);
        vec mina_outer   = titanlib::util::subset(values_mina, indices_global_outer);
        vec minv_outer   = titanlib::util::subset(values_minv, indices_global_outer);
        vec maxa_outer   = titanlib::util::subset(values_maxa, indices_global_outer);
        vec maxv_outer   = titanlib::util::subset(values_maxv, indices_global_outer);


        // -~- Compute the background 
        vec bvalues_outer = background( background_elab_type, elevs_outer, values_outer, num_min_prof, min_elev_diff, value_minp, value_maxp, background_values, indices_global_outer, debug);

        if(debug) std::cout << "... background ok ..." << std::endl;

        // -~- If deviations between backgrounds and observations are small then flag = 0

        /* if observations and background are almost identical (within the range of valid estimates),
           then take a shortcut and flag = 0 for all observations in the inner circle */

        int j = indices_outer_test[0];
        if(bvalues_outer[j] > minv_outer[j] && bvalues_outer[j] < maxv_outer[j]) { 
            int g = indices_global_test[0];
            flags[g] = 0;
            if (debug) std::cout << " small_innov - index " << g << std::endl;
            continue;
        }

        // -~- SCT on a selection of observations to test in the inner circle

        res = sct_core( lats_outer, lons_outer, elevs_outer, values_outer, bvalues_outer, min_horizontal_scale, max_horizontal_scale, kth_closest_obs_horizontal_scale, vertical_scale, eps2_outer, value_minp, value_maxp, mina_outer, maxa_outer, minv_outer, maxv_outer, tpos_outer, tneg_outer, indices_global_test, indices_outer_inner, indices_outer_test, indices_inner_test, debug, na, true, thrown_out, scores, flags);

        // problems during the matrix inversion
        if ( !res) {
            flags[curr] = 100;
            if(debug) std::cout << " oi - flags=100 - index " << curr << std::endl; 
            continue;
        }

        count_oi++; // one more OI

    }  // end loop over observations

    std::cout << "QC missing - Removing " << thrown_out << " observations. Number of OI " << count_oi << std::endl;
    double e_time0 = titanlib::util::clock();
    std::cout << e_time0 - s_time0 << " secs" << std::endl;

    /* final check on the bad observations
       it may happen that a good observation is flagged as a bad one because of the order
       the SCT has been done and the uncertainty made in the estimation of the background 
       (remember that bad observations may have been used to get the background).
       This final check is made only on bad observations and uses only good observations. */

    if(debug) std::cout << " +++++ final check on the bad observations ++++++++++++++++" << std::endl;
    s_time0 = titanlib::util::clock();
        
    // reset this number each loop (this is for breaking if we don't throw anything new out)
    thrown_out = 0; 
        
    // diagnostic. count the number of times OI is performed
    count_oi = 0;

    // loop over observations
    for(int curr=0; curr < p; curr++) {

        if(debug) std::cout << "===> curr " << curr << " ===============" << std::endl;

        // -~- Can the observation be used as a centroid for the definition of inner and outer circles?
        if(obs_test[curr] != 1 || flags[curr] != 1 || values[curr] < value_minp || values[curr] > value_maxp) {
            if(debug) std::cout << "..skip " << curr << std::endl;
            continue;
        }

        // -~- Define outer/inner circles and which are the observations to test

        // get all neighbours that are close enough (inside outer circle)
        vec distances;
        ivec indices_global_outer_guess = tree.get_neighbours_with_distance(lats[curr], lons[curr], outer_radius, num_max_outer, true, distances);

        // set all the indices linking the different levels

        ivec indices_global_outer, indices_global_test, indices_outer_inner, indices_outer_test, indices_inner_test;
        bool res = set_indices( indices_global_outer_guess, obs_test, flags, distances, inner_radius, curr, indices_global_outer, indices_global_test, indices_outer_inner, indices_outer_test, indices_inner_test);
            
        int p_outer = indices_global_outer.size();
        int p_inner = indices_outer_inner.size();
        int p_test = indices_global_test.size();
        if(debug) std::cout << "p_outer inner test " << p_outer << " " << p_inner << " " << p_test << std::endl;
            
        // -~- Decide if there are enough observations in the outer/inner circles

        if(p_outer < num_min_outer) {
            flags[curr] = 12;
            if(debug) std::cout << "@@isolated (outer) " << curr << std::endl;
            continue; 
        }

        if( p_inner < 2) {
            flags[curr] = 11;
            if(debug) std::cout << "@@isolated (inner) " << curr << std::endl;
            continue;
        }
            
        // -~- Set vectors on the outer circle

        vec lons_outer   = titanlib::util::subset(lons, indices_global_outer);
        vec elevs_outer  = titanlib::util::subset(elevs, indices_global_outer);
        vec lats_outer   = titanlib::util::subset(lats, indices_global_outer);
        vec values_outer = titanlib::util::subset(values, indices_global_outer);
        vec eps2_outer   = titanlib::util::subset(eps2, indices_global_outer);
        vec tpos_outer   = titanlib::util::subset(tpos, indices_global_outer);
        vec tneg_outer   = titanlib::util::subset(tneg, indices_global_outer);
        vec mina_outer   = titanlib::util::subset(values_mina, indices_global_outer);
        vec minv_outer   = titanlib::util::subset(values_minv, indices_global_outer);
        vec maxa_outer   = titanlib::util::subset(values_maxa, indices_global_outer);
        vec maxv_outer   = titanlib::util::subset(values_maxv, indices_global_outer);

        // Debug for indices
        if(debug) {
            ivec flags_outer;
            flags_outer.resize( p_outer, na);
            for(int i=0; i < p_outer; i++) {
               int g = indices_global_outer[i];
               flags_outer[i] = flags[g];
            }
            std::cout << "curr lats lons obs:" << std::endl;
            std::cout << curr << " " << lats[curr] << " " << lons[curr] << " " << values[curr] << std::endl;
            std::cout << "indices_global_outer - i g lats lons obs:" << std::endl;
            for(int i=0; i<p_outer; i++) {
                int g = indices_global_outer[i];
                std::cout << std::setprecision(6) << i << " " << g << " " << lats[g] << " " << lons[g] << " " << values[g] << " " << flags[g] << std::endl;
            }
            std::cout << "outer - i lats lons obs:" << std::endl;
            for(int i=0; i<p_outer; i++) {
                std::cout << std::setprecision(6) << i << " " << lats_outer[i] << " " << lons_outer[i] << " " << values_outer[i] << " " << flags_outer[i] << std::endl;
            }
            std::cout << "indices_outer_inner - l i lats lons obs:" << std::endl;
            for(int l=0; l<p_inner; l++) {
                int i = indices_outer_inner[l];
                std::cout << std::setprecision(6) << l << " " << i << " " << lats_outer[i] << " " << lons_outer[i] << " " << values_outer[i] << " " << flags_outer[i] << std::endl;
            }
            std::cout << "indices_outer_test - l m lats lons obs:" << std::endl;
            for(int m=0; m<p_test; m++) {
                int i = indices_outer_test[m];
                std::cout << std::setprecision(6) << m << " " << i << " " << lats_outer[i] << " " << lons_outer[i] << " " << values_outer[i] << " " << flags_outer[i] << std::endl;
            }
            std::cout << "indices_inner_test - l m lats lons obs:" << std::endl;
            for(int m=0; m<p_test; m++) {
                int l = indices_inner_test[m];
                int i = indices_outer_inner[l];
                std::cout << std::setprecision(6) << m << " " << l << " " << lats_outer[i] << " " << lons_outer[i] << " " << values_outer[i] << " " << flags_outer[i] << std::endl;
            }
        }

        // -~- Compute the background 
        vec bvalues_outer = background( background_elab_type, elevs_outer, values_outer, num_min_prof, min_elev_diff, value_minp, value_maxp, background_values, indices_global_outer, debug);

        if(debug) { 
            int j = indices_outer_test[0];
            for(int i=0; i < p_outer; i++) {
                std::cout << std::setprecision(6) << "values backg minv maxv " << values_outer[i] << " " << bvalues_outer[i] << " " << minv_outer[i] << " " << maxv_outer[i] << " " << std::endl; 
            }
            std::cout << "j " << j << std::endl; 
        }
        if(debug) std::cout << "... background ok ..." << std::endl;

        // -~- If deviations between backgrounds and observations are small then flag = 0

        /* if observations and background are almost identical (within the range of valid estimates),
           then take a shortcut and flag = 0 for all observations in the inner circle */

        int j = indices_outer_test[0];
        if(bvalues_outer[j] > minv_outer[j] && bvalues_outer[j] < maxv_outer[j]) { 
            int g = indices_global_test[0];
            flags[g] = 0;
            if (debug) std::cout << " small_innov - index " << g << std::endl;
            continue;
        }

        // -~- SCT within the inner circle
        // note that the only observation to test is curr
        res = sct_core( lats_outer, lons_outer, elevs_outer, values_outer, bvalues_outer, min_horizontal_scale, max_horizontal_scale, kth_closest_obs_horizontal_scale, vertical_scale, eps2_outer, value_minp, value_maxp, mina_outer, maxa_outer, minv_outer, maxv_outer, tpos_outer, tneg_outer, indices_global_test, indices_outer_inner, indices_outer_test, indices_inner_test, debug, na, true, thrown_out, scores, flags);

        // problems during the matrix inversion
        if ( !res) {
            flags[curr] = 100;
            if(debug) std::cout << " oi - flags=100 - index " << curr << std::endl; 
            continue;
        }

        count_oi++; // one more OI

    }  // end loop over observations

    std::cout << "Re-check bad obs - Removing " << thrown_out << " observations. Number of OI " << count_oi << std::endl;
    e_time0 = titanlib::util::clock();
    std::cout << e_time0 - s_time0 << " secs" << std::endl;

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
            } else if( flags[curr] == 11) {
                count_iso_inner++;
            } else if( flags[curr] == 12) {
                count_iso_outer++;
            } else if( flags[curr] == 100) {
                count_fail_matinv++;
            } else {
                count_impossible++;
            }
        }
        std::cout << std::setprecision(3) << "summary - # TOT good bad missing isolated(inner) isolated(outer): " << p << " " << count_good << " " << count_bad << " " << count_missing << " " << count_iso_inner << " " << count_iso_outer << std::endl; 
        if(count_fail_matinv > 0) 
            std::cout << std::setprecision(3) << "!!!! failure in matrix inversion: " << count_fail_matinv << std::endl; 
        if(count_impossible > 0) 
            std::cout << std::setprecision(3) << "!!!! unknown flag: " << count_impossible << std::endl; 
    }
    //
    
    std::cout << ">> Total Time " << e_time0 - s_time << "secs" << std::endl;

    //
    return flags;
}
// end SCT //


//=============================================================================

/*----------------------------------------------------------------------------
 BACKGROUND FUNCTIONS */

//+ Compute the background for all the observations in the outer circle
vec background( std::string background_elab_type,
                const vec& elevs, 
                const vec& values, 
                int num_min_prof, 
                float min_elev_diff,
                float value_minp,
                float value_maxp,
                const vec& external_background_values,
                const ivec& indices_global_outer,
                bool debug) {
    vec bvalues;
    int p = values.size();
    
    if( background_elab_type == "vertical_profile") {
        // vertical profile is independent from the curr-th observation
        bvalues = compute_vertical_profile(elevs, elevs, values, num_min_prof, min_elev_diff, debug);
    } else if( background_elab_type == "vertical_profile_Theil_Sen") {
        // vertical profile is independent from the curr-th observation
        bvalues = compute_vertical_profile_Theil_Sen(elevs, elevs, values, num_min_prof, min_elev_diff, debug);
    } else if( background_elab_type == "mean_outer_circle"){
        // vertical profile is independent from the curr-th observation
        float mean_val = std::accumulate(values.begin(), values.end(), 0.0) / p;
        bvalues.resize(p, mean_val);
    } else if( background_elab_type == "median_outer_circle"){
        // vertical profile is independent from the curr-th observation
        float median_val = titanlib::util::compute_quantile( 0.5, values);
        bvalues.resize(p, median_val);
    } else if( background_elab_type == "external"){
        bvalues = titanlib::util::subset(external_background_values, indices_global_outer);
    }
    // set the range 
    for(int i=0; i<p; i++) {
        if(bvalues[i] < value_minp) bvalues[i] = value_minp;
        if(bvalues[i] > value_maxp) bvalues[i] = value_maxp;
    }
    //
    return bvalues;
}

vec compute_vertical_profile_Theil_Sen(const vec& elevs, const vec& oelevs, const vec& values, int num_min_prof, double min_elev_diff, bool debug) {
    // Starting value guesses
    double gamma = -0.0065;
    double min_diff = (double) min_elev_diff;
    double meanT = std::accumulate(values.begin(), values.end(), 0.0) / values.size();

    // Special case when all observations have the same elevation
    if ( std::min_element(elevs.begin(),elevs.end()) == std::max_element(elevs.begin(),elevs.end())) {
        vec vp( oelevs.size(), meanT);
        return vp;
    }

    // Check if terrain is too flat
    double z05 = titanlib::util::compute_quantile(0.05, elevs);
    double z95 = titanlib::util::compute_quantile(0.95, elevs);

    // should we use the basic or more complicated vertical profile?
    bool use_basic = elevs.size() < num_min_prof || (z95 - z05) < min_diff;

    // Theil-Sen (Median-slope) Regression (Wilks (2019), p. 284)
    int n = elevs.size();
    int no = oelevs.size();
    int nm = n * (n - 1) / 2;
    std::vector<float> q(n), m(nm);
    double m_median;
    if ( use_basic) {
        m_median = gamma;
    } else {
        int k = 0;
        for(int i=0; i<(n-1); i++) {
            for(int j=(i+1); j<n; j++) {
                float dz = fabs( elevs[i] - elevs[j]);
                if ( dz < 1 ) {
                    m[k] = 0;
                } else {
                    m[k] = ( values[i] - values[j]) / ( elevs[i] - elevs[j]);
                }
                k++;
            }
        }
        m_median = titanlib::util::compute_quantile( 0.5, m);
    }
    for(int i=0; i<n; i++) {
        q[i] = values[i] - m_median * elevs[i];
    }
    double q_median = titanlib::util::compute_quantile( 0.5, q);
    vec vp(no);
    if (debug) std::cout << "Theil-Sen vp - m q " << m_median  << " " << q_median << std::endl;
    for(int i=0; i<no; i++) {
      vp[i] = q_median + m_median * oelevs[i];
    }
    //
    return vp;
}

vec compute_vertical_profile(const vec& elevs, const vec& oelevs, const vec& values, int num_min_prof, double min_elev_diff, bool debug) {
    // Starting value guesses
    double gamma = -0.0065;
    double min_diff = (double) min_elev_diff;
    double a = 5.0;

    double meanT = std::accumulate(values.begin(), values.end(), 0.0) / values.size();

    // Special case when all observations have the same elevation
    if ( std::min_element(elevs.begin(),elevs.end()) == std::max_element(elevs.begin(),elevs.end())) {
        vec vp( oelevs.size(), meanT);
        return vp;
    }
    double exact_p10 = titanlib::util::compute_quantile(0.10, elevs);
    double exact_p90 = titanlib::util::compute_quantile(0.90, elevs);
    
    // optimize inputs for VP (using Nelder-Mead Simplex algorithm)
    const gsl_multimin_fminimizer_type *T = gsl_multimin_fminimizer_nmsimplex2;
    gsl_multimin_fminimizer *s = NULL;
    gsl_vector *ss;
    gsl_multimin_function vp_optim;


    int iter = 0;
    int status;
    double size;

    // data (params) that needs to be passed into vp
    double nd = (double) elevs.size(); // cast + resize to double
    double nod = (double) oelevs.size(); // cast + resize to double
    std::vector<double> delevs(elevs.begin(), elevs.end());
    std::vector<double> doelevs(oelevs.begin(), oelevs.end());
    std::vector<double> dvalues(values.begin(), values.end());
    double * dpelevs = delevs.data();
    double * dpoelevs = doelevs.data();
    double * dpvalues = dvalues.data();
    double * data[3] = {&nd, dpelevs, dpvalues};

    // Check if terrain is too flat
    double z05 = titanlib::util::compute_quantile(0.05, elevs);
    double z95 = titanlib::util::compute_quantile(0.95, elevs);

    // should we use the basic or more complicated vertical profile?
    bool use_basic = elevs.size() < num_min_prof || (z95 - z05) < min_diff;

    gsl_vector* input;
    if(use_basic) {
        vp_optim.n = 2;
        input = gsl_vector_alloc(vp_optim.n);
        gsl_vector_set(input, 0, meanT);
        gsl_vector_set(input, 1, gamma);
        vp_optim.f = basic_vertical_profile_optimizer_function;
    }
    else {
        vp_optim.n = 5;
        input = gsl_vector_alloc(vp_optim.n);
        gsl_vector_set(input, 0, meanT);
        gsl_vector_set(input, 1, gamma);
        gsl_vector_set(input, 2, a);
        gsl_vector_set(input, 3, exact_p10);
        gsl_vector_set(input, 4, exact_p90);
        vp_optim.f = vertical_profile_optimizer_function;
    }
    ss = gsl_vector_alloc (vp_optim.n);

    gsl_vector_set_all (ss, 1.0);
    vp_optim.params = data;

    s = gsl_multimin_fminimizer_alloc (T, vp_optim.n);
    gsl_multimin_fminimizer_set (s, &vp_optim, input, ss);
    do {
        iter++;
        status = gsl_multimin_fminimizer_iterate(s);

        if (status)
            break;

        size = gsl_multimin_fminimizer_size (s);
        status = gsl_multimin_test_size (size, 1e-2);
    }
    while (status == GSL_CONTINUE && iter < 100);

    // then actually calculate the vertical profile using the minima
    vec vp;
    if(use_basic) {
        vp = basic_vertical_profile(nod, dpoelevs, gsl_vector_get(s->x, 0), gsl_vector_get(s->x, 1));
        if(debug) std::cout << "vp - meanT=" << gsl_vector_get(s->x, 0) << " gamma=" << gsl_vector_get(s->x, 1) << std::endl;
    }
    else {
        // then actually calculate the vertical profile using the minima
        vp = vertical_profile(nod, dpoelevs,  gsl_vector_get(s->x, 0), gsl_vector_get(s->x, 1),
                gsl_vector_get(s->x, 2), gsl_vector_get(s->x, 3), gsl_vector_get(s->x, 4));
        if(debug) std::cout << "vp - t0 " << gsl_vector_get(s->x, 0) << " gamma=" << gsl_vector_get(s->x, 1) << " a=" << gsl_vector_get(s->x, 2) << " h0=" << gsl_vector_get(s->x, 3) << " h1i=" << gsl_vector_get(s->x, 4) << std::endl;
    }

    gsl_vector_free(input);
    gsl_vector_free(ss);
    gsl_multimin_fminimizer_free (s);

    return vp;
}

double basic_vertical_profile_optimizer_function(const gsl_vector *v, void *data) {
    double **p = (double **)data;
    int n = (int) *p[0]; // is of type double but should be an int
    double *dpelevs = p[1];
    double *dpvalues = p[2];

    // the parameters to mess with
    double meanT = gsl_vector_get(v,0);
    double gamma = gsl_vector_get(v,1);

    // give everything to vp to compute t_out
    vec t_out = basic_vertical_profile(n, dpelevs, meanT, gamma);
    // RMS
    float total = 0;
    for(int i=0; i<n; i++) {
        total += pow((t_out[i]-dpvalues[i]),2);
    }
    double value = log(pow((total / n),0.5));
    return value;
}

vec basic_vertical_profile(const int n, const double *elevs, const double t0, const double gamma) {
    vec t_out(n, -999);
    for(int i=0; i<n; i++)
        t_out[i] = t0 + gamma*elevs[i];
    return t_out;
}

double vertical_profile_optimizer_function(const gsl_vector *v, void *data) {
    double **p = (double **)data;
    int n = (int) *p[0]; // is of type double but should be an int
    double *dpelevs = p[1];
    double *dpvalues = p[2];

    // the parameters to mess with
    double meanT = gsl_vector_get(v,0);
    double gamma = gsl_vector_get(v,1);
    double a = gsl_vector_get(v,2);
    double exact_p10 = gsl_vector_get(v,3);
    double exact_p90 = gsl_vector_get(v,4);

    // give everything to vp to compute t_out
    vec t_out = vertical_profile(n, dpelevs, meanT, gamma, a, exact_p10, exact_p90);
    // RMS
    double total = 0;
    for(int i=0; i<n; i++) {
        total += pow((t_out[i]-dpvalues[i]),2);
    }
    double value = log(pow((total / n),0.5));
    return value;
}

vec vertical_profile(const int n, const double *elevs, const double t0, const double gamma, const double a, const double h0, const double h1i) {
    double h1 = h0 + fabs(h1i); // h1<-h0+abs(h1i)
    // loop over the array of elevations
    vec t_out;
    t_out.resize(n, -999);

    for(int i=0; i<n; i++) {
        // define some bools
        bool z_le_h0 = elevs[i] <= h0; // z.le.h0<-which(z<=h0)
        bool z_ge_h1 = elevs[i] >= h1; // z.ge.h1<-which(z>=h1)
        bool z_in = (elevs[i]>h0 && elevs[i]<h1); // z.in<-which(z>h0 & z<h1)
        if(z_le_h0) {
            t_out[i] = t0+gamma*elevs[i]-a;
        }
        if(z_ge_h1) {
            t_out[i] = t0+gamma*elevs[i];
        }
        if(z_in) {
            t_out[i] = t0+gamma*elevs[i]-a/2*(1+cos(M_PI*(elevs[i]-h0)/(h1-h0)));
        }
    }
    return t_out;
}


//----------------------------------------------------------------------------//
// HELPER FUNCTIONS

//+ Matrix inversion based on LU-decomposition
bool invert_matrix(const boost::numeric::ublas::matrix<float>& input, boost::numeric::ublas::matrix<float>& inverse) {
    using namespace boost::numeric::ublas;
    typedef permutation_matrix<std::size_t> pmatrix;
    // create a working copy of the input
    matrix<float> A(input);
    // create a permutation matrix for the LU-factorization
    pmatrix pm(A.size1());
    // perform LU-factorization
    int res = lu_factorize(A,pm);
    if( res != 0 ) return false;
    // create identity matrix of "inverse"
    inverse.assign(boost::numeric::ublas::identity_matrix<float>(A.size1()));
    // backsubstitute to get the inverse
    lu_substitute(A, pm, inverse);
    return true;
}

//+ Set indices
bool set_indices( const ivec& indices_global_outer_guess, 
                  const ivec& obs_test,
                  const ivec& flags, 
                  const vec& dist_outer_guess, 
                  float inner_radius, 
                  int test_just_this,
                  ivec& indices_global_outer, 
                  ivec& indices_global_test, 
                  ivec& indices_outer_inner, 
                  ivec& indices_outer_test, 
                  ivec& indices_inner_test) {
/*-----------------------------------------------------------------------------
 Set all the vectors of indices linking the vectors operating over the
 different regions that we have defined: global, outer circle, inner circle,
 observations to test within the inner circle
 If the index test_just_this is set to a value greater than 0 (i.e. it 
 does point to an observation of the global vectors), then the only
 observation to test is test_just_this.
 
 Rules for an observation to be assigned to:
 1) outer circle. Within outer_radius from the centroid and flag!=1
 2) inner circle. In the outer circle AND within inner_radius from the centroid
 3) test. In the inner circle AND flag!=0 AND the user want to test it 

----------------------------------------------------------------------------*/
    int p_outer_guess = indices_global_outer_guess.size();
    indices_global_outer.reserve(p_outer_guess);
    indices_global_test.reserve(p_outer_guess);
    indices_outer_inner.reserve(p_outer_guess);
    indices_outer_test.reserve(p_outer_guess);
    indices_inner_test.reserve(p_outer_guess);
    int i = -1; // counter for indices on outer
    int l = -1; // counter for indices on inner
    for(int count=0; count<p_outer_guess; count++) {
        int g = indices_global_outer_guess[count];
        if(flags[g] != 1 && g != test_just_this) {
            indices_global_outer.push_back(g);
            i++;
            if(dist_outer_guess[count] <= inner_radius) {
                indices_outer_inner.push_back(i);
                l++;
                if (test_just_this < 0 && flags[g] != 0 && obs_test[g] == 1) {
                    indices_global_test.push_back(g);
                    indices_outer_test.push_back(i);
                    indices_inner_test.push_back(l);
                }
            }
        }
    }
    if ( test_just_this >= 0) {
        // add global index to global_outer indices as its last element
        indices_global_outer.push_back(test_just_this);
        // add global index to global_test as its only element
        indices_global_test.push_back(test_just_this);
        // adjust outer_inner, last element of outer_inner points to last element of outer 
        indices_outer_inner.push_back(indices_global_outer.size()-1);
        // adjust outer_test, only element of test points to last element of outer 
        indices_outer_test.push_back(indices_global_outer.size()-1);
        // adjust inner_test, only element of test points to last element of inner 
        indices_inner_test.push_back(indices_outer_inner.size()-1);
    }
    //
    return true;
}

//----------------------------------------------------------------------------//
// SCT WITHIN THE INNER CIRCLE

// Smallets SCT unit
bool sct_core( const vec& lats, 
               const vec& lons, 
               const vec& elevs, 
               const vec& yo, 
               const vec& yb, 
               float Dh_min, 
               float Dh_max, 
               int kth_close, 
               float Dz, 
               const vec& eps2, 
               float minp,
               float maxp,
               const vec& mina, 
               const vec& maxa, 
               const vec& minv, 
               const vec& maxv, 
               const vec& tpos, 
               const vec& tneg, 
               const ivec& indices_global_test, 
               const ivec& indices_outer_inner, 
               const ivec& indices_outer_test, 
               const ivec& indices_inner_test, 
               bool debug, 
               float na, 
               bool set_flag0, 
               int& thrown_out,
               vec& scores,
               ivec& flags) {
/*-----------------------------------------------------------------------------
 Spatial consistency test for a selection of observations to test within the
 inner circle, considering all the observations in the outer circle.

 The inner and outer circles are defined with respect to a centroid observation.

 Each observation to test is compared against an estimated value (cvanalysis)
 obtained considering the neighbouring observations but without considering the
 observation under test (leave-one-out cross-validation). In addition, the SCT
 evaluate the likelihood that the observed value is draw from a Gaussian(good)  
 rather than a uniform(bad) PDF. The evaluation requires the computation of an
 additional estimate, the analysis, which is the best estimate of the value
 at an observation location when we consider all the observed values, even 
 the observation under test. Both analysis and cvanalysis are calculated
 taking into consideration all observations in the outer circle.
 
 SCT scores:
 chi = sqrt( analysis_residual * cvanalysis_residual) 
   analysis residual   = yo - ya
   cvanalysis residual = yo - yav
 z = (chi - mu) / (sigma + sigma_mu)
   mu = median( chi(*))
   sigma = inter-quartile range of chi(*)
   sigma_mu = sigma / sqrt(n(*))
 (*) statistics based on the n observations that are within the inner circle
 and with a cvanalysis within the range of admissible values.

 SCT fails (the test shows the uniform PDF as the most likely) when:
 z > tpos, yo >= yav
 z > tneg, yo <  yav

 Our prior knowledge includes the definition of two confidence levels for the 
 cvanalysis:
  - valid values, range [minv,maxv]
      cvanalysis is so close that it validates the corresponding observation
  - admissible values, range [mina,maxa] 
      if cvanalysis is outside this range, then the observation is bad 
  Note that the statistics defining the spread used in the SCT is based on 
  deviations between observations and estimates that lie within the range of 
  admissible values only.

 Decision tree:

 ALL observations are set to bad when:
  the cvanalyses all lie outside the range of admissible values 

 ONE observation is deemed as bad when:
  fails the SCT with the largest z-score 
  AND 
  cvanalysis is outside the range of valid values

 ALL observations are set to good when:
  the user allows the procedure to flag them as good ones (set_flag0 is true)
  AND
   the cvanalyses all lie within the range of valid values
   OR
   none of the cvanalysis fail the SCT test

  Output values:

   flags, global vector. return updated quality control flag values (1 = bad, 0 = good)
   
   scores, global vector. return updated scores. Note that only for the case of 
                          a single bad observation we store the z-score
   
   thrown_out, integer. update the number of observations flagged as bad ones 

-----------------------------------------------------------------------------*/

    using namespace boost::numeric::ublas;

    // init
    int p_outer = yo.size();
    int p_inner = indices_outer_inner.size();
    int p_test = indices_global_test.size();

    /* Compute Dh. 
     The location-dependent horizontal de-correlation lenght scale 
     used for the background error correlation matrix */
    boost::numeric::ublas::matrix<float> disth(p_outer, p_outer);
    boost::numeric::ublas::matrix<float> distz(p_outer, p_outer);
    boost::numeric::ublas::vector<float> Dh(p_outer);

    for(int i=0; i < p_outer; i++) {
        vec Dh_vector(p_outer);
        int k = 0;
        for(int j=i; j < p_outer; j++) {
            disth(i, j) = titanlib::util::calc_distance( lats[i], lons[i], lats[j], lons[j]);
            distz(i, j) = fabs( elevs[i] - elevs[j]);
            if(i != j) {
                disth(j, i) = disth(i, j);
                distz(j, i) = distz(i, j);
                Dh_vector[k] = disth(i, j);
                k++;
            }
        }
        for(int j=0; j <= i; j++) {
            Dh_vector[k] = disth(i, j);
            k++;
        }
        // find distance to the k-th closest observations
        Dh(i) = titanlib::util::findKclosest( kth_close, Dh_vector); 
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
    
    // Compute S + eps2*I and store it in S 
    boost::numeric::ublas::matrix<float> S(p_outer,p_outer);
    boost::numeric::ublas::matrix<float> Sinv(p_outer,p_outer);
    for(int i=0; i < p_outer; i++) {
        for(int j=i; j < p_outer; j++) {
            float value = std::exp(-.5 * std::pow((disth(i, j) / Dh_mean), 2) - .5 * std::pow((distz(i, j) / Dz), 2));
            if(i==j) { // weight the diagonal, this also ensure an invertible matrix
                value = value + eps2[i];
            } else {
                S(j,i) = value;
            }
            S(i,j) = value;
        }
    }
    
    // Compute ( S + eps2 * I )^(-1) ans store it in Sinv
    bool b = invert_matrix(S, Sinv);
    if( !b)  return( false);

    // Definitions
    boost::numeric::ublas::vector<float> Sinv_d(p_outer);

    // Matrix multiplications ( S + eps2 * I )^(-1) * (yo - yb)
    for(int i=0; i<p_outer; i++) {
        S(i,i) -= eps2[i];
        float acc = 0;
        for(int j=0; j<p_outer; j++) {
            acc += Sinv(i,j) * (yo[j] - yb[j]);
        }
        Sinv_d(i) = acc;
    }

    // obtain chi = sqrt( analysis_residual * cvanalysis_residual) 
    // two chi vectors are obtained:
    // 1. chi_inner, chi for all observations in the inner circle
    // 2. chi_stat, chi values used to extract summary stastics 
    //    observations in the inner circle AND yav has an admissible value   
    vec chi_stat;
    vec chi_stat_alt;
    chi_stat.reserve(p_inner);
    chi_stat_alt.reserve(p_inner);
    boost::numeric::ublas::vector<float> chi_inner(p_inner);
    boost::numeric::ublas::vector<float> chi_inner_alt(p_inner);
    boost::numeric::ublas::vector<float> yav(p_inner);
    if(debug) std::cout << "analysis - i l lon lat yo yb ya yav chi chi_alt" << std::endl;
    for(int l=0; l<p_inner; l++) {
        int i = indices_outer_inner[l];
        float ya = yb[i];
        for(int j=0; j<p_outer; j++) 
            ya += S(i,j)*Sinv_d(j);
        yav(l) = yo[i] - Sinv_d(i) / Sinv(i,i);
        if(ya < minp) ya = minp;
        if(ya > maxp) ya = maxp;
        if(yav(l) < minp) yav(l) = minp;
        if(yav(l) > maxp) yav(l) = maxp;
        chi_inner(l) = std::sqrt( (yo[i] - ya) * (yo[i] - yav(l))); 
        chi_inner_alt(l) = std::sqrt( eps2[i] / ( 1 + eps2[i])) * ( maxv[i] - minv[i]);
        if ( yav(l) >= mina[i] && yav(l) <= maxa[i]) { 
            chi_stat.push_back(chi_inner(l));
            chi_stat_alt.push_back(chi_inner_alt(l));
        }
        if(debug) std::cout << std::setprecision(6) << i << " " << l << " " << lons[i] << " " << lats[i] << " " << yo[i] << " " << yb[i] << " " << ya << " " << yav(l) << " " << chi_inner(l) << " " << chi_inner_alt(l) << std::endl;
    }
    
    // chi_stat is empty (yav all outside the range of admissible values), set all flags to 1
    if(chi_stat.size() == 0) {
        for(int m=0; m<p_test; m++) {
            int g = indices_global_test[m];
            flags[g] = 1;
            thrown_out++;
            if(debug) std::cout << std::setprecision(3) << "chi_stat empty - flag as bad " << g << std::endl;
        }
        return true;       
    }

    // chi summary statistics
    float mu = titanlib::util::compute_quantile( 0.5, chi_stat);
    float sigma = titanlib::util::compute_quantile( 0.75, chi_stat) - titanlib::util::compute_quantile( 0.25, chi_stat);
    float sigma_alt = titanlib::util::compute_quantile( 0.75, chi_stat_alt) - titanlib::util::compute_quantile( 0.25, chi_stat_alt);
    if (sigma_alt > sigma) sigma = sigma_alt;
 
    if(sigma == 0) {
        if(debug) std::cout << "sigma = 0, not possible to perform the spatial check " << std::endl;
        return true;       
    }
    float sigma_mu = sigma / std::sqrt( chi_stat.size());

    if(debug) {
        std::cout << " chi_stat, dim = " << chi_stat.size() << std::endl;
        for(int l=0; l<chi_stat.size(); l++) 
            std::cout << std::setprecision(6) << chi_stat[l] << std::endl;
        std::cout << std::setprecision(6) << " mu sigma sigma_mu " << mu << " " << sigma << " " << sigma_mu << std::endl;
    }
   
    // z = ( chi - mu ) / ( sigma + sigma_mu )
    float zmx = -10000;
    int mmx = -1;
    if(debug) std::cout << "z_test - i l lons lats yo yb yav chi z " << std::endl;
    for(int m=0; m<p_test; m++) {
        int i = indices_outer_test[m];
        int l = indices_inner_test[m];
        float z = (chi_inner(l) - mu) / (sigma + sigma_mu);
        if (z > zmx && (yav(l) < minv[i] || yav(l) > maxv[i])) {
          zmx = z;
          mmx = m;
        } 
        if(debug) std::cout << std::setprecision(6) << i << " " << l << " " << lons[i] << " " << lats[i] << " " << yo[i] << " " << yb[i] << " " << yav(l) << " " << chi_inner(l) << " " << z << std::endl;
    }

    // Decision making 
    int i = indices_outer_test[mmx];
    int l = indices_inner_test[mmx];
    // set the threshold
    float thr = tneg[i];
    if ( (yo[i] - yav(l)) >= 0) thr = tpos[i];
    // if the largest z is larger than the threshold, then flag IT as bad ...
    if ( zmx > thr) {
        int g = indices_global_test[mmx];
        scores[g] = zmx;
        flags[g] = 1;
        thrown_out++;
        if(debug) std::cout << std::setprecision(6) << "SCT failed - flag as bad - index yo yav chi z " << g << " " << yo[i] << " " << yav(l) << " " << chi_inner(l) << " " << zmx << std::endl;
    // ... BUT if the largest z is smaller than the threshold (and we are confindet in flagging), then flag ALL as good
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
