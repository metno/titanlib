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
// for the diagnostics file
#include <iostream>
#include <fstream>  // For file handling
#include <cstdio>      // For std::remove

using namespace titanlib;

// helpers

// Remove flagged elements from a vector of indices
void sct_dual_with_fg_remove_flagged(ivec& indices, vec& distances, const ivec& flags1, const ivec& flags2, const int& index_to_keep, const int& nmax, const int& f1, const int& f2);

// Leave-One-Out Integral Data Influece based on Optimal Interpolation
void sct_dual_with_fg_cvidi( vec& cvidis_yes, vec& cvidis_no, float& Dh_mean, const vec& lons, const vec& lats, const vec& elevs, const vec& vyes, const vec& byes, const vec& eps2, const float& Dz, const float& Dh_min, const float& Dh_max );

//-----------------------------------------------------------------------------
// start SCT with first guess //
ivec titanlib::sct_dual_with_fg(const Points& points,
        const vec& values,
        const vec& background_values,
        const vec& event_thresholds,
        ConditionType condition,
        int num_min,
        int num_max,
        float inner_radius,
        float outer_radius,
        int num_iterations,
        float min_horizontal_scale,
        float max_horizontal_scale,
        float vertical_scale,
        const vec& eps2,
        bool diagnostics,
        const std::string& filename_diagnostics,
        vec& sct_cvidis_yes,
        vec& sct_cvidis_no,
        const ivec& obs_to_check) {
//-----------------------------------------------------------------------------
    // Diagnostics
    std::ofstream outfile; 
    if (diagnostics) {
        // Check if the file exists and delete it
        if (std::ifstream(filename_diagnostics)) {
            std::remove(filename_diagnostics.c_str());
            std::cout << "File existed and was deleted.\n";
        }
        outfile.open(filename_diagnostics, std::ios_base::app); // Open file in append mode
        if (!outfile) { // Check if the file opened successfully
            std::cerr << "Error opening file!" << std::endl;
        }
        outfile << "it;loop;curr;i;index;lon;lat;z;vyes;byes;dh;flags_d;cvidi_yes_d;cvidi_no_d;flags_r;cvidi_yes_r;cvidi_no_r;saved_r;flags;sct_cvidi_yes;sct_cvidi_no;" << std::endl;
    }

    // Constants
    // code for "Not Available". Any type of missing data
    float na = -999.; 
    // Min cvidis required to flag, avoiding overly isolated observations
    float min_cvidi = 0.1; 
    // Min required cvidis difference to flag, ensuring conservative flagging and accounting for uncertainty
    float cvidi_buffer = 0.05; 

    // Initializations
    const vec& lats = points.get_lats();
    const vec& lons = points.get_lons();
    const vec& elevs = points.get_elevs();
    const int s = values.size();
    ivec flags(s, 0);
    vec r(s, na);
    sct_cvidis_yes.clear();
    sct_cvidis_yes.resize(s, na);
    sct_cvidis_no.clear();
    sct_cvidis_no.resize(s, na);

    // Safe checks
    if(lats.size() != s)
        throw std::runtime_error("Lats does not have same size as values");
    if(lons.size() != s)
        throw std::runtime_error("Lons does not have same size as values");
    if(elevs.size() != s)
        throw std::runtime_error("Elevs does not have same size as values");
    if(eps2.size() != s)
        throw std::runtime_error("Eps2 does not have same size as values");
    if(background_values.size() != s)
        throw std::runtime_error("Background_values does not have same size as values");
    if(obs_to_check.size() > 0 && obs_to_check.size() != s)
        throw std::invalid_argument("obs_to_check must empty or have the same size as values");
    if(num_min < 2)
        throw std::invalid_argument("num_min must be > 1");
    if(num_max < num_min)
        throw std::invalid_argument("num_max must be > num_min");
    if(num_iterations < 1)
        throw std::invalid_argument("num_iterations must be >= 1");
    if(min_horizontal_scale <= 0)
        throw std::invalid_argument("min_horizontal_scale must be > 0");
    if(max_horizontal_scale < min_horizontal_scale)
        throw std::invalid_argument("max_horizontal_scale must be >= min_horizontal_scale");
    if(vertical_scale <= 0)
        throw std::invalid_argument("vertical_scale must be > 0");
    if(inner_radius < 0)
        throw std::invalid_argument("inner_radius must be >= 0");
    if(outer_radius < inner_radius)
        throw std::invalid_argument("outer_radius must be >= inner_radius");
    for(int i = 0; i < eps2.size(); i++) {
        if(eps2[i] <= 0)
            throw std::invalid_argument("All eps2 values must be > 0");
    }

    // k-d trees help for nearest neighbor searches
    titanlib::KDTree tree(lats, lons);

    // Flag stations without elevation
    for(int curr=0; curr < s; curr++) {
        if(!titanlib::is_valid(elevs[curr])) {
            flags[curr] = 1;
        }
    }

    // Set event thresholds
    if (event_thresholds.size() == s) {
        std::copy(event_thresholds.begin(), event_thresholds.end(), r.begin());
    } else {
        std::fill(r.begin(), r.end(), event_thresholds[0]);
    }

    // Populate vectors with event occurrences (1 = yes, 0 = no)
    vec vyes(s, 0);
    vec byes(s, 0);
    for(int curr=0; curr < s; curr++) {
        if ( condition == titanlib::Eq) {
            if ( values[curr] == r[curr]) vyes[curr] = 1;
            if ( background_values[curr] == r[curr]) byes[curr] = 1;
        } else if ( condition == titanlib::Geq) {
            if ( values[curr] >= r[curr]) vyes[curr] = 1;
            if ( background_values[curr] >= r[curr]) byes[curr] = 1;
        } else if ( condition == titanlib::Gt) {
            if ( values[curr] >  r[curr]) vyes[curr] = 1;
            if ( background_values[curr] > r[curr]) byes[curr] = 1;
        } else if ( condition == titanlib::Leq) {
            if ( values[curr] <= r[curr]) vyes[curr] = 1;
            if ( background_values[curr] <= r[curr]) byes[curr] = 1;
        } else if ( condition == titanlib::Lt) {
            if ( values[curr] <  r[curr]) vyes[curr] = 1;
            if ( background_values[curr] < r[curr]) byes[curr] = 1;
        }
    }
    
    //---------------------------------------------------------------------
    // SCT dual ITERATION
    for(int iteration = 0; iteration < num_iterations; iteration++) {

        double s_time0 = titanlib::clock();

        if(diagnostics) {
            std::cout << "==========================================" << std::endl;
            std::cout << ">> ITERATION: " << iteration << std::endl;
        }

        int thrown_out = 0; // reset this counter in each loop; used to break if no new elements are removed
        ivec checked_d(s, 0);  // Tracks which observations have been checked
        ivec flags_d(s, 0);  // Flags from Detection Loop
        ivec flags_r(s, 0);  // Flags from Stray data redemption loop

        // Flagging becomes progressively more challenging with each iteration 
        // to prevent excessive removal along the borders between "yes" and "no" regions
        float cvidi_buffer_it = (iteration+1) * cvidi_buffer; 

        // diagnostics
        vec cvidis_yes_d(s, na); // Scores from Detection Loop
        vec cvidis_no_d(s, na); // Scores from Detection Loop
        vec cvidis_yes_r(s, na); // Scores from Stray data redemption loop
        vec cvidis_no_r(s, na); // Scores from Stray data redemption loop
        vec dh_d(s, na); // Dh from Detection Loop 
        vec dh_r(s, na); // Dh from Stray data redemption loop
        vec saved_r(s, 0); // Flagged by Detection Loop & Cluster Preservation loop, saved by Stray data redemption loop

        //---------------------------------------------------------------------
        // DETECTION LOOP
        // Identifying suspect observations
        
        // Diagnostics: Count the number of calls to the OI function 
        // (this is the most time-consuming part)
        int count_oi_d = 0;
        double s_time = titanlib::clock();
        for(int curr=0; curr < s; curr++) {
            
            // break out if observation already checked
            if(checked_d[curr] > 0) 
                continue;
            
            // break out if observation not to check OR already flagged
            if((obs_to_check[curr] != 1) || (flags[curr] != 0)) {
                checked_d[curr] = 1;
                continue;
            }

            // get all neighbours that are close enough and not flagged in previous iterations
            vec distances;
            ivec neighbour_indices = tree.get_neighbours_with_distance(lats[curr], lons[curr], outer_radius, distances, true);
            sct_dual_with_fg_remove_flagged(neighbour_indices, distances, flags, flags, curr, num_max, 0, 0);
            int s_box = neighbour_indices.size();
            
            // break out if observation isolated (note that we don't flag it)
            if(s_box < num_min) {
                checked_d[curr] = 1;
                continue; 
            }

            // call SCT dual with this box 
            vec lons_box = titanlib::subset(lons, neighbour_indices);
            vec elevs_box = titanlib::subset(elevs, neighbour_indices);
            vec lats_box = titanlib::subset(lats, neighbour_indices);
            vec vyes_box = titanlib::subset(vyes, neighbour_indices);
            vec byes_box = titanlib::subset(byes, neighbour_indices);
            vec eps2_box = titanlib::subset(eps2, neighbour_indices);

            // perform statistical interpolation within the box
            vec cvidis_yes_box(s_box, na);
            vec cvidis_no_box(s_box, na);
            float Dh_mean = 0.;
            sct_dual_with_fg_cvidi( cvidis_yes_box, cvidis_no_box, Dh_mean, lons_box, lats_box, elevs_box, vyes_box, byes_box, eps2_box, vertical_scale, min_horizontal_scale, max_horizontal_scale);
            count_oi_d++;

            // save variables for diagnostics
            dh_d[curr] = Dh_mean;
            
            // detect suspect observations within the inner circle
            // NOTE: In the detection loop, an observation can be checked multiple times
            //  because it may fell into several inner circles
            //  if it is flagged as suspect just once, then it is not checked again (in this loop)
            for(int i = 0; i < s_box; i++) {
                int index = neighbour_indices[i];
                if(obs_to_check[index] != 1) {
                    checked_d[index] = 1;
                    continue;
                }
                if( (distances[i] <= inner_radius) && 
                    ( (cvidis_yes_box[i] >= min_cvidi) || 
                      (cvidis_no_box[i] >= min_cvidi)) ){ 
                    // condition to identify a candidate for flagging
                    if( ((vyes_box[i] == 0) && (cvidis_yes_box[i] > (cvidis_no_box[i] + cvidi_buffer_it))) || 
                        ((vyes_box[i] == 1) && (cvidis_no_box[i]  > (cvidis_yes_box[i] + cvidi_buffer_it))) ) {
                        flags_d[index] = 1;
                        cvidis_yes_d[index] = cvidis_yes_box[i];
                        cvidis_no_d[index] = cvidis_no_box[i];
                    } else {
                        // in case of multiple checking of the same observation, then keep the worst cvidis
                        if ((vyes_box[i] == 0) && (cvidis_no_box[i] < cvidis_no_d[index])) {
                            cvidis_yes_d[index] = cvidis_yes_box[i];
                            cvidis_no_d[index] = cvidis_no_box[i];
                        } else if ((vyes_box[i] == 1) && (cvidis_yes_box[i] < cvidis_yes_d[index])) {
                            cvidis_yes_d[index] = cvidis_yes_box[i];
                            cvidis_no_d[index] = cvidis_no_box[i];
                        }
                    }
                    checked_d[index] = 1;
                }
            }

            // diagnostics
            if(diagnostics) {
                for(int i = 0; i < s_box; i++) {
                    int index = neighbour_indices[i];
                    // "it;loop;curr;i;index;lon;lat;z;vyes;byes;dh;flags_d;cvidi_yes_d;cvidi_no_d;flags_r;cvidi_yes_r;cvidi_no_r;saved_r;flags;sct_cvidi_yes_r;sct_cvidi_no_r;"
                    outfile << iteration << ";1;" << curr << ";" << i << ";" << index << ";";
                    outfile << std::fixed << std::setprecision(5) << lons_box[i] << ";" << lats_box[i] << ";";
                    outfile << std::fixed << std::setprecision(0) << elevs_box[i] << ";";
                    outfile << std::fixed << std::setprecision(2) << vyes_box[i] << ";" << byes_box[i] << ";";
                    outfile << std::fixed << std::setprecision(0) << dh_d[curr] << ";";
                    outfile << std::fixed << std::setprecision(0) << flags_d[index] << ";";
                    outfile << std::fixed << std::setprecision(2) << cvidis_yes_d[index] << ";";
                    outfile << std::fixed << std::setprecision(2) << cvidis_no_d[index] << ";";
                    outfile << std::fixed << std::setprecision(0) << "-999;-999;-999;-999;-999;-999;-999;" << std::endl; 
                }
            }

        } // End DETECTION LOOP

        // diagnostics
        if(diagnostics) {
            double e_time = titanlib::clock();
            std::cout << "Detection Loop, Time " << std::fixed << std::setprecision(0) << e_time - s_time << " sec" << std::endl;
        }

        //---------------------------------------------------------------------
        // STRAY DATA REDEMPTION LOOP
        // Bringing back good observations that got caught in the wrong crowd

        // Diagnostics: Count the number of calls to the OI function 
        // (this is the most time-consuming part)
        int count_oi_r = 0;
        s_time = titanlib::clock();
        for(int curr=0; curr < s; curr++) {

            // break out if observation has NOT been flagged in the cluster preservation loop
            if(flags_d[curr] == 0) 
                continue;

            // get all neighbours that are close enough (consider only non-flagged neighbours)
            vec distances;
            ivec neighbour_indices = tree.get_neighbours_with_distance(lats[curr], lons[curr], outer_radius, distances, true);
            sct_dual_with_fg_remove_flagged(neighbour_indices, distances, flags, flags_d, curr, num_max, 0, 0);
            int s_box = neighbour_indices.size();

            // break out if observation isolated
            // NOTE: this time we flag it (different from detection loop) 
            //       because it has been flagged by the previous loops
            if(s_box < num_min) {
                flags_r[curr] = 1;
                continue; 
            }

            // identify the index of curr observation in the box 
            int curr_in_box = 0;
            for(int i = 0; i < s_box; i++) {
                if(neighbour_indices[i] == curr) { 
                    curr_in_box = i;
                    break;
                }
            }

            // call SCT dual with this box 
            vec lons_box = titanlib::subset(lons, neighbour_indices);
            vec elevs_box = titanlib::subset(elevs, neighbour_indices);
            vec lats_box = titanlib::subset(lats, neighbour_indices);
            vec vyes_box = titanlib::subset(vyes, neighbour_indices);
            vec byes_box = titanlib::subset(byes, neighbour_indices);
            vec eps2_box = titanlib::subset(eps2, neighbour_indices);

            // perform statistical interpolation within the box
            vec cvidis_yes_box(s_box, na);
            vec cvidis_no_box(s_box, na);
            float Dh_mean = 0.;
            sct_dual_with_fg_cvidi( cvidis_yes_box, cvidis_no_box, Dh_mean, lons_box, lats_box, elevs_box, vyes_box, byes_box, eps2_box, vertical_scale, min_horizontal_scale, max_horizontal_scale);
            count_oi_r++;

            // save variables for diagnostics
            dh_r[curr] = Dh_mean;

            // check the curr observation
            if( ( (cvidis_yes_box[curr_in_box] >= min_cvidi) || (cvidis_no_box[curr_in_box] >= min_cvidi)) && 
                ( ((vyes_box[curr_in_box] == 0) && (cvidis_yes_box[curr_in_box] > (cvidis_no_box[curr_in_box] + cvidi_buffer_it))) || 
                  ((vyes_box[curr_in_box] == 1) && (cvidis_no_box[curr_in_box]  > (cvidis_yes_box[curr_in_box] + cvidi_buffer_it)))) ) {
                flags_r[curr] = 1;
            } else {
                saved_r[curr] = 1;
            }
            cvidis_yes_r[curr] = cvidis_yes_box[curr_in_box];
            cvidis_no_r[curr] = cvidis_no_box[curr_in_box];

            if(diagnostics) {
                for(int i = 0; i < s_box; i++) {
                    int index = neighbour_indices[i];
                    // "it;loop;curr;i;index;lon;lat;z;vyes;byes;dh;flags_d;cvidi_yes_d;cvidi_no_d;flags_r;cvidi_yes_r;cvidi_no_r;saved_r;flags;sct_cvidi_yes_r;sct_cvidi_no_r;"
                    outfile << iteration << ";2;" << curr << ";" << i << ";" << index << ";";
                    outfile << std::fixed << std::setprecision(5) << lons_box[i] << ";" << lats_box[i] << ";";
                    outfile << std::fixed << std::setprecision(0) << elevs_box[i] << ";";
                    outfile << std::fixed << std::setprecision(2) << vyes_box[i] << ";" << byes_box[i] << ";";
                    outfile << std::fixed << std::setprecision(0) << dh_d[curr] << ";";
                    outfile << std::fixed << std::setprecision(0) << flags_d[index] << ";";
                    outfile << std::fixed << std::setprecision(2) << cvidis_yes_d[index] << ";";
                    outfile << std::fixed << std::setprecision(2) << cvidis_no_d[index] << ";";
                    outfile << std::fixed << std::setprecision(0) << flags_r[index] << ";";
                    outfile << std::fixed << std::setprecision(2) << cvidis_yes_r[index] << ";";
                    outfile << std::fixed << std::setprecision(2) << cvidis_no_r[index] << ";";
                    outfile << std::fixed << std::setprecision(0) << saved_r[index] << ";";
                    outfile << std::fixed << std::setprecision(0) << "-999;-999;-999;" << std::endl; 
                }
            }

        } // End STRAY DATA REDEMPTION LOOP

        // diagnostics
        if(diagnostics) {
            double e_time = titanlib::clock();
            std::cout << "Stray Data Redemption Loop, Time " << std::fixed << std::setprecision(0) << e_time - s_time << " sec" << std::endl;
        }

        //---------------------------------------------------------------------
        // FLAG ASSIGNMENT STEP
        // Definitive flag is assigned to each observation based on the outcomes 
        // of the detection, cluster preservation, and redemption processes
        int tot_checked_d = 0;
        int tot_flagged_d = 0;
        int tot_flagged_r = 0;
        int tot_saved_r = 0;
        int tot_thrown_out = 0;
        for(int curr=0; curr < s; curr++) {
            // Flagged is a previous SCT iteration
            if(flags[curr] != 0) {
                tot_thrown_out++;
                if(diagnostics && checked_d[curr] == 1)
                    tot_checked_d++;
                continue;
            }
            // Fresh new flagged observation
            if(flags_r[curr] == 1) {
                flags[curr] = 1;
                // Keep the last cvidis 
                sct_cvidis_yes[curr] = cvidis_yes_r[curr];
                sct_cvidis_no[curr] = cvidis_no_r[curr];
                thrown_out++;
                tot_thrown_out++;
            // Not flagged observation
            } else if (saved_r[curr] == 1) {
                sct_cvidis_yes[curr] = cvidis_yes_r[curr];
                sct_cvidis_no[curr] = cvidis_no_r[curr];
            } else {
                sct_cvidis_yes[curr] = cvidis_yes_d[curr];
                sct_cvidis_no[curr] = cvidis_no_d[curr];
            }
            // diagnostics
            if(diagnostics) {
                if(checked_d[curr] == 1)
                    tot_checked_d++;
                if(flags_d[curr] == 1)
                    tot_flagged_d++;
                if(flags_r[curr] == 1)
                    tot_flagged_r++;
                if(saved_r[curr] == 1)
                    tot_saved_r++;
            }
        } // End FLAG ASSIGNMENT STEP

        if(diagnostics) {
            for(int curr=0; curr < s; curr++) {
                // "it;loop;curr;i;index;lon;lat;z;vyes;byes;dh;flags_d;cvidi_yes_d;cvidi_no_d;flags_r;cvidi_yes_r;cvidi_no_r;saved_r;flags;sct_cvidi_yes_r;sct_cvidi_no_r;"
                outfile << iteration << ";3;" << curr << ";0;" << curr << ";";
                outfile << std::fixed << std::setprecision(5) << lons[curr] << ";" << lats[curr] << ";";
                outfile << std::fixed << std::setprecision(0) << elevs[curr] << ";";
                outfile << std::fixed << std::setprecision(2) << vyes[curr] << ";" << byes[curr] << ";";
                outfile << std::fixed << std::setprecision(0) << "-999;";
                outfile << std::fixed << std::setprecision(0) << flags_d[curr] << ";";
                outfile << std::fixed << std::setprecision(2) << cvidis_yes_d[curr] << ";";
                outfile << std::fixed << std::setprecision(2) << cvidis_no_d[curr] << ";";
                outfile << std::fixed << std::setprecision(0) << flags_r[curr] << ";";
                outfile << std::fixed << std::setprecision(2) << cvidis_yes_r[curr] << ";";
                outfile << std::fixed << std::setprecision(2) << cvidis_no_r[curr] << ";";
                outfile << std::fixed << std::setprecision(0) << saved_r[curr] << ";";
                outfile << std::fixed << std::setprecision(0) << flags[curr] << ";";
                outfile << std::fixed << std::setprecision(2) << sct_cvidis_yes[curr] << ";";
                outfile << std::fixed << std::setprecision(2) << sct_cvidis_no[curr] << ";" << std::endl; 
            }
            std::cout << "-- Detection Loop -- " << std::endl;
            std::cout << "# obs " << s << " Checked " << tot_checked_d << " Flagged " << tot_flagged_d << " Number of OI loop " << count_oi_d << std::endl;
            std::cout << "-- Stray Data Redemption Loop -- " << std::endl;
            std::cout << "# checked " << tot_flagged_d << " Flagged " << tot_flagged_r << " Saved " << tot_saved_r << " Number of OI loop " << count_oi_r << std::endl;
            std::cout << "-- Summary -- " << std::endl;
            std::cout << "Removing " << thrown_out << " (Total removed " << tot_thrown_out << " , " << std::setprecision(0) << 100*tot_thrown_out/s << "%)" << std::endl;
            double e_time0 = titanlib::clock();
            std::cout << "Time " << std::fixed << std::setprecision(0) << e_time0 - s_time0 << " sec" << std::endl;
        }

        //---------------------------------------------------------------------
        // Check if we can STOP
        if(thrown_out == 0) {
            if(iteration + 1 < num_iterations) {
                if(diagnostics)
                    std::cout << "Stopping early after " << iteration + 1<< " iterations" << std::endl;
            }
            break;
        }

    } // End loop over iterations

    return flags;
}
// end SCT with FG //

//----------------------------------------------------------------------------//
// HELPER FUNCTIONS

// Remove flagged elements from a vector of indices
void sct_dual_with_fg_remove_flagged(
    ivec& indices, 
    vec& distances, 
    const ivec& flags1, 
    const ivec& flags2, 
    const int& index_to_keep, 
    const int& nmax, 
    const int& f1, 
    const int& f2
) {
//-----------------------------------------------------------------------------
    ivec indices_new;
    vec distances_new;

    indices_new.reserve(indices.size());
    distances_new.reserve(indices.size());

    for(int i=0; i<indices.size(); i++) {
        if((flags1[indices[i]] == f1 && flags2[indices[i]] == f2) || indices[i] == index_to_keep) {
            indices_new.push_back(indices[i]);
            distances_new.push_back(distances[i]);
        }
    }

    if(indices_new.size() > nmax) {
        int N = indices_new.size();
        std::vector<std::pair<float,int> > pairs(N);
        for(int i = 0; i < indices_new.size(); i++) {
            pairs[i] = std::pair<float, int>(distances_new[i], indices_new[i]);
        }
        std::sort(pairs.begin(), pairs.end(), titanlib::sort_pair_first<float,int>());
        distances_new.clear();
        indices_new.clear();
        distances_new.resize(nmax);
        indices_new.resize(nmax);
        for(int i = 0; i < nmax; i++) {
            distances_new[i] = pairs[i].first;
            indices_new[i] = pairs[i].second;
        }
    }

    indices = indices_new;
    distances = distances_new;

} // END sct_dual_with_fg_remove_flagged

// Leave-One-Out Integral Data Influece based on Optimal Interpolation
void sct_dual_with_fg_cvidi(
    vec& cvidis_yes, 
    vec& cvidis_no, 
    float& Dh_mean, 
    const vec& lons, 
    const vec& lats, 
    const vec& elevs, 
    const vec& vyes, 
    const vec& byes, 
    const vec& eps2, 
    const float& Dz, 
    const float& Dh_min, 
    const float& Dh_max 
) {
//-----------------------------------------------------------------------------
    int s_dim = vyes.size();

    // Categorize observations into "yes" and "no" groups based on vyes values
    ivec indices_yes, indices_no;
    indices_yes.reserve(s_dim);
    indices_no.reserve(s_dim);
    for (int i = 0; i < s_dim; i++) {
        (vyes[i] == 1 ? indices_yes : indices_no).push_back(i);
    }
    int sy_dim = indices_yes.size(), sn_dim = indices_no.size();

   // If all observations belong to one category, set flags accordingly
    if (sy_dim == 0 || sn_dim == 0) {
        std::fill(cvidis_yes.begin(), cvidis_yes.end(), sy_dim > 0);
        std::fill(cvidis_no.begin(), cvidis_no.end(), sn_dim > 0);
        return;
    }

    // Compute distances
    boost::numeric::ublas::matrix<float> disth(s_dim, s_dim, 0.0f), distz(s_dim, s_dim, 0.0f);
    for (int i = 0; i < s_dim; i++) {
        for (int j = i + 1; j < s_dim; j++) {
            double dist_h = titanlib::calc_distance(lats[i], lons[i], lats[j], lons[j]);
            double dist_z = fabs(elevs[i] - elevs[j]);
            disth(i, j) = disth(j, i) = dist_h;
            distz(i, j) = distz(j, i) = dist_z;
        }
    }

    // Compute horizontal decorrelation length
    if (Dh_min == Dh_max) {
        Dh_mean = Dh_min;
    } else {
        vec Dh_vector(s_dim - 1);
        boost::numeric::ublas::vector<float> Dh(s_dim);
        for (int i = 0; i < s_dim; i++) {
            int k = 0;
            for (int j = 0; j < s_dim; j++) {
                if (i != j) Dh_vector[k++] = disth(i, j);
            }
            Dh(i) = titanlib::compute_quantile(0.10, Dh_vector);
        }
        Dh_mean = std::accumulate(Dh.begin(), Dh.end(), 0.0f) / Dh.size();
        Dh_mean = std::max(Dh_min, std::min(Dh_mean, Dh_max));
    }

    // Construct background error correlation matrix S
    boost::numeric::ublas::matrix<float> S(s_dim, s_dim);
    double inv_Dh_mean = 1.0 / Dh_mean;
    double inv_Dz = 1.0 / Dz;
    double neg_half = -0.5;
    for(int i=0; i < s_dim; i++) {
        for(int j=i; j < s_dim; j++) {
            double norm_dist_h = disth(i, j) * inv_Dh_mean;
            double norm_dist_z = distz(i, j) * inv_Dz;
            double corr_value = std::exp(neg_half * (norm_dist_h * norm_dist_h + norm_dist_z * norm_dist_z));
            S(i, j) = corr_value;
            S(j, i) = corr_value; // Use symmetry
        }
    }

    // Loop for cvidis_no (k=0) and cvidis_yes (k=1)
    for (int k = 0; k < 2; k++) {
        const ivec &indices_k = (k == 0 ? indices_no : indices_yes);
        const ivec &indices_l = (k == 0 ? indices_yes : indices_no);
        int sk_dim = (k == 0 ? sn_dim : sy_dim), sl_dim = (k == 0 ? sy_dim : sn_dim);

        // Extract submatrices
        boost::numeric::ublas::matrix<float> Gl(sl_dim, sk_dim), Sk(sk_dim, sk_dim);
        for (int i = 0; i < sk_dim; ++i) {
            for (int j = 0; j < sk_dim; ++j) {
                Sk(i, j) = S(indices_k[i], indices_k[j]);
            }
            for (int j = 0; j < sl_dim; ++j) {
                Gl(j, i) = S(j, indices_k[i]);
            }
        }

        // Initialize innovation vector dk 
        boost::numeric::ublas::vector<float> dk(sk_dim, 1.0f);

        // Compute Sk = S+R
        for (int i = 0; i < sk_dim; i++) {
            int j = indices_k[i];
            Sk(i, i) += (vyes[j] != byes[j]) ? 1.0f : eps2[j];
        }

        // Invert (S + R) matrix
        boost::numeric::ublas::matrix<float> Skinv(sk_dim, sk_dim);
        titanlib::invert_matrix(Sk, Skinv);

        // (S+R)^(-1) * (yo-yb)
        boost::numeric::ublas::vector<float> Skinv_d(sk_dim);
        Skinv_d = boost::numeric::ublas::prod(Skinv, dk);

        // Define cvidis
        boost::numeric::ublas::vector<float> cvidis(s_dim);
        
        // Compute cvidis for l-locations (where we don't have observations)
        boost::numeric::ublas::vector<float> cvidis_l = boost::numeric::ublas::prod(Gl, Skinv_d);
        for (int i = 0; i < sl_dim; i++) {
            int j = indices_l[i];
            cvidis(j) = cvidis_l(i);
        }

        // Compute cvidis for k-locations (where we have observations)
        if (sk_dim == 1) {
            for (int i = 0; i < sk_dim; ++i) {
                cvidis(indices_k[i]) = 0;
            }
        } else {
            for (int i = 0; i < sk_dim; ++i) {
                cvidis(indices_k[i]) = 1 - Skinv_d(i) / Skinv(i, i);
            }
        }

        // Store results in cvidis_yes or cvidis_no
        auto &target = (k == 0) ? cvidis_no : cvidis_yes;
        std::copy(cvidis.begin(), cvidis.end(), target.begin());
    }

} // END sct_dual_with_fg_oi
