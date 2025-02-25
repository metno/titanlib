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
#include <boost/math/distributions/students_t.hpp>
// for the diagnostics file
#include <iostream>
#include <fstream>  // For file handling
//#include <filesystem>  // C++17 for checking file existence and deleting
#include <cstdio>      // For std::remove

using namespace titanlib;

// helpers
void sct_with_fg_remove_flagged(ivec& indices, vec& distances, const ivec& flags1, const ivec& flags2, const int& index_to_keep, const int& nmax, const int& f1, const int& f2);
void sct_with_fg_oi(vec& ares, vec& cvres, vec& innov, float& Dh_mean, const vec& lons, const vec& lats, const vec& elevs, const vec& values, const vec& backg, const vec& eps2, const float& Dz, const float& Dh_min, const float& Dh_max, const float& vmin, const float& vmax);
void sct_with_fg_oi_multi(vec& ares, vec& cvres, vec& innov, float& Dh_mean, const vec& lons, const vec& lats, const vec& elevs, const vec& values, const vec& backg, const vec& eps2, const float& Dz, const float& Dh_min, const float& Dh_max, const float& vmin, const float& vmax, const int& n_multi);
float sct_with_fg_sig2_estimate(const vec& ares, const vec& innov, const float& vmin, const float& q);
float sct_with_fg_sig2_estimate_alt(const vec& values, const float& vmin);

// start SCT with first guess //
ivec titanlib::sct_with_fg(const Points& points,
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
        float max_horizontal_scale,
        float vertical_scale,
        const vec& pos,
        const vec& neg,
        const vec& eps2,
        const vec& min_obs_var,
        bool diagnostics,
        const std::string& filename_diagnostics,
        vec& sct_score,
        const ivec& obs_to_check) {

    //diagnostics
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
        outfile << "it;loop;curr;i;index;lon;lat;z;yo;yb;ya;yav;dh;sig2;tc;flags_d;score_d;flags_c;score_c;saved_c;flags_r;score_r;saved_r;flags;sct_score;" << std::endl;
    }

//    float alpha = 0.05; // Significance level (5%) used in the Cluster Preservation Loop
    float alpha = 0.01; // Significance level (5%) used in the Cluster Preservation Loop

    const vec& lats = points.get_lats();
    const vec& lons = points.get_lons();
    const vec& elevs = points.get_elevs();

    const int s = values.size();
    if(lats.size() != s)
        throw std::runtime_error("Lats does not have same size as values");
    if(lons.size() != s)
        throw std::runtime_error("Lons does not have same size as values");
    if(elevs.size() != s)
        throw std::runtime_error("Elevs does not have same size as values");
    if(pos.size() != s)
        throw std::runtime_error("Pos does not have same size as values");
    if(neg.size() != s || eps2.size() != s)
        throw std::runtime_error("Neg does not have same size as values");
    if(eps2.size() != s)
        throw std::runtime_error("Eps2 does not have same size as values");
    if(min_obs_var.size() != s)
        throw std::runtime_error("min_obs_var does not have same size as values");
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
    if(vertical_scale <= 0)
        throw std::invalid_argument("vertical_scale must be > 0");
    if(inner_radius < 0)
        throw std::invalid_argument("inner_radius must be >= 0");
    if(outer_radius < inner_radius)
        throw std::invalid_argument("outer_radius must be >= inner_radius");
    for(int i = 0; i < eps2.size(); i++) {
        if(eps2[i] <= 0)
            throw std::invalid_argument("All eps2 values must be > 0");
        if(min_obs_var[i] <= 0)
            throw std::invalid_argument("All min_obs_var values must be > 0");
        if(pos[i] < 0)
            throw std::invalid_argument("All pos values must be >= 0");
        if(neg[i] < 0)
            throw std::invalid_argument("All neg values must be >= 0");
    }

    ivec flags(s, 0);
    sct_score.clear();
    sct_score.resize(s, 0);

    titanlib::KDTree tree(lats, lons);

    // Flag stations without elevation
    for(int curr=0; curr < s; curr++) {
        if(!titanlib::is_valid(elevs[curr])) {
            flags[curr] = 1;
        }
    }
    
    // Screen observations to check before SCT iteration
    ivec obs_to_check_before_sct(s, 1);
    for(int curr=0; curr < s; curr++) {
        if(obs_to_check[curr] != 1) {
            obs_to_check_before_sct[curr] = obs_to_check[curr];
            continue;
        }
        // break out if observation close to first guess
        float buffer = 2 * std::sqrt(min_obs_var[curr]);
        if((values[curr] >= (background_values[curr]-buffer)) && (values[curr] <= (background_values[curr]+buffer))) {
            obs_to_check_before_sct[curr] = 0;
        }
    }

    //---------------------------------------------------------------------
    // SCT ITERATION
    for(int iteration = 0; iteration < num_iterations; iteration++) {
        double s_time0 = titanlib::clock();

        int thrown_out = 0; // reset this number each loop (this is for breaking if we don't throw anything new out)
        int saved = 0; // reset this number each loop (this is for knowing how many observations have been given a 2nd chance)

        ivec checked_d(s, 0);  // Keep track of which observations have been checked
        ivec flags_d(s, 0);  // 
        ivec flags_c(s, 0);  // 
        ivec flags_r(s, 0);  // 
        vec score_d(s, 0); // 
        vec score_c(s, 0); // 
        vec score_r(s, 0); // 
        // diagnostics
        vec dh_d(s, 0); 
        vec dh_c(s, 0); 
        vec dh_r(s, 0); 
        vec sig2_d(s, 0); 
        vec sig2_c(s, 0); 
        vec sig2_r(s, 0); 
        vec saved_c(s, 0); 
        vec saved_r(s, 0); 

        //---------------------------------------------------------------------
        // DETECTION LOOP
        // Identifying suspect observations
        int count_oi_d = 0;
        for(int curr=0; curr < s; curr++) {
            
            // break out if observation already checked
            if(checked_d[curr] > 0) 
                continue;
            
            // break out if observation not ot check OR already flagged
            if((obs_to_check_before_sct[curr] != 1) || (flags[curr] != 0)) {
                checked_d[curr] = 1;
                continue;
            }

            // get all neighbours that are close enough
            vec distances;
            ivec neighbour_indices = tree.get_neighbours_with_distance(lats[curr], lons[curr], outer_radius, distances, true);
            sct_with_fg_remove_flagged(neighbour_indices, distances, flags, flags, curr, num_max, 0, 0);
            int s_box = neighbour_indices.size();
            
            // break out if observation isolated (note that we don't flag it)
            if(s_box < num_min) {
                checked_d[curr] = 1;
                continue; 
            }

            // call SCT with this box 
            vec lons_box = titanlib::subset(lons, neighbour_indices);
            vec elevs_box = titanlib::subset(elevs, neighbour_indices);
            vec lats_box = titanlib::subset(lats, neighbour_indices);
            vec values_box = titanlib::subset(values, neighbour_indices);
            vec backg_box = titanlib::subset(background_values, neighbour_indices);
            vec eps2_box = titanlib::subset(eps2, neighbour_indices);

            // perform statistical interpolation within the box
            vec ares(s_box,0);
            vec cvres(s_box,0);
            vec innov(s_box,0);
            float Dh_mean = 0.;
            sct_with_fg_oi_multi(ares, cvres, innov, Dh_mean, lons_box, lats_box, elevs_box, values_box, backg_box, eps2_box, vertical_scale, min_horizontal_scale, max_horizontal_scale, values_min, values_max, 1);
            count_oi_d++;
            // estimate observation error variance
            float sig2o = sct_with_fg_sig2_estimate(ares, innov, min_obs_var[curr], 0.5); 

            // save variables for diagnostics
            dh_d[curr] = Dh_mean;
            sig2_d[curr] = sig2o;
            
            // detect suspect observations within the inner circle
            // NOTE: In the detection loop, an observation can be checked multiple times
            //  because it may fell into several inner circles
            //  if it is flagged as suspect just once, then it is not checked again (in this loop)
            for(int i = 0; i < s_box; i++) {
                int index = neighbour_indices[i];
                if(obs_to_check_before_sct[index] != 1) {
                    checked_d[index] = 1;
                    continue;
                }
                float dist = distances[i];
                if(dist <= inner_radius) {
                    float sig2o_i = std::max( sig2o, min_obs_var[index]);
                    float score = cvres[i] * ares[i] / sig2o_i;
                    assert(titanlib::is_valid(score));
                    // condition to identify a candidate for flagging
                    if((cvres[i] < 0 && score > pos[index]) || (cvres[i] >= 0 && score > neg[index])) {
                        sig2_d[index] = sig2o_i;
                        flags_d[index] = 1;
                    }
                    // in case of multiple checking of the same observation, then keep the highest score
                    // (used in the cluster preservation loop)
                    score_d[index] = std::max( score_d[index], score);
                    checked_d[index] = 1;
                }
            }

            // diagnostics
            if(diagnostics) {
                for(int i = 0; i < s_box; i++) {
                    int index = neighbour_indices[i];
                    // "it;loop;curr;i;index;lon;lat;z;yo;yb;ya;yav;dh;sig2;flags_d;score_d;flags_c;score_c;saved_c;flags_r;score_r;saved_r;flags;sct_score;"
                    outfile << iteration << ";1;" << curr << ";" << i << ";" << index << ";";
                    outfile << std::fixed << std::setprecision(5) << lons_box[i] << ";" << lats_box[i] << ";";
                    outfile << std::fixed << std::setprecision(0) << elevs_box[i] << ";";
                    outfile << std::fixed << std::setprecision(2) << values_box[i] << ";" << backg_box[i] << ";";
                    outfile << std::fixed << std::setprecision(2) << ares[i] + values_box[i] << ";" << cvres[i] + values_box[i] << ";";
                    outfile << std::fixed << std::setprecision(0) << dh_d[curr] << ";";
                    outfile << std::fixed << std::setprecision(3) << sig2_d[curr] << ";";
                    outfile << std::fixed << std::setprecision(0) << "-999;";
                    outfile << std::fixed << std::setprecision(0) << flags_d[index] << ";";
                    outfile << std::fixed << std::setprecision(3) << score_d[index] << ";";
                    outfile << std::fixed << std::setprecision(0) << "-999;-999;-999;-999;-999;-999;-999;-999;" << std::endl; 
                }
            }

        } // End DETECTION LOOP

        //---------------------------------------------------------------------
        // CLUSTER PRESERVATION LOOP
        // Saving data that blends well with neighbors
        int count_oi_c = 0;
        for(int curr=0; curr < s; curr++) {

            // break out if observation has NOT been flagged in the detection loop
            if(flags_d[curr] == 0) 
                continue;

            // get all neighbours that are close enough (exact same neighbours as in the detection loop)
            vec distances;
            ivec neighbour_indices = tree.get_neighbours_with_distance(lats[curr], lons[curr], outer_radius, distances, true);
            sct_with_fg_remove_flagged(neighbour_indices, distances, flags, flags, curr, num_max, 0, 0);
            int s_box = neighbour_indices.size();
            // break out if observation isolated (note that we flag it)
//            if(s_box < num_min) {
//                flags_c[curr] = 1;
//                continue; 
//            }

            // identify the index of curr observation in the box 
            int curr_in_box = 0;
            for(int i = 0; i < s_box; i++) {
                int index = neighbour_indices[i];
                if(index == curr) { 
                    curr_in_box = i;
                    break;
                }
            }

            // call SCT with this box 
            vec lons_box = titanlib::subset(lons, neighbour_indices);
            vec elevs_box = titanlib::subset(elevs, neighbour_indices);
            vec lats_box = titanlib::subset(lats, neighbour_indices);
            vec values_box = titanlib::subset(values, neighbour_indices);
            vec backg_box = titanlib::subset(background_values, neighbour_indices);
            vec eps2_box = titanlib::subset(eps2, neighbour_indices);
            
            // perform statistical interpolation within the box
            vec ares(s_box,0);
            vec cvres(s_box,0);
            vec innov(s_box,0);
            float Dh_mean = 0.;
            sct_with_fg_oi_multi(ares, cvres, innov, Dh_mean, lons_box, lats_box, elevs_box, values_box, backg_box, eps2_box, vertical_scale, min_horizontal_scale, max_horizontal_scale, values_min, values_max, 4);
            count_oi_c++;
            // estimate variance in the scores within the box
            float sig2o = sig2_d[curr]; 

            // save variables for diagnostics
            dh_c[curr] = Dh_mean;
            sig2_c[curr] = sig2o;

            // re-check
            float score = cvres[curr_in_box] * ares[curr_in_box] / sig2o;
            assert(titanlib::is_valid(score));
            score_c[curr] = score;
            if((cvres[curr_in_box] < 0 && score > pos[curr]) || (cvres[curr_in_box] >= 0 && score > neg[curr])) {
                flags_c[curr] = 1;
            } else {
                saved_c[curr] = 1;
            }

            // diagnostics
            if(diagnostics) {
                for(int i = 0; i < s_box; i++) {
                    int index = neighbour_indices[i];
                    // "it;loop;curr;i;index;lon;lat;z;yo;yb;ya;yav;dh;sig2;flags_d;score_d;flags_c;score_c;saved_c;flags_r;score_r;saved_r;flags;sct_score;"
                    outfile << iteration << ";2;" << curr << ";" << i << ";" << index << ";";
                    outfile << std::fixed << std::setprecision(5) << lons_box[i] << ";" << lats_box[i] << ";";
                    outfile << std::fixed << std::setprecision(0) << elevs_box[i] << ";";
                    outfile << std::fixed << std::setprecision(2) << values_box[i] << ";" << backg_box[i] << ";";
                    outfile << std::fixed << std::setprecision(2) << ares[i] + values_box[i] << ";" << cvres[i] + values_box[i] << ";";
                    outfile << std::fixed << std::setprecision(0) << dh_c[curr] << ";";
                    outfile << std::fixed << std::setprecision(3) << sig2_c[curr] << ";";
                    outfile << std::fixed << std::setprecision(3) << "-999;"; 
                    outfile << std::fixed << std::setprecision(0) << flags_d[index] << ";";
                    outfile << std::fixed << std::setprecision(3) << score_d[index] << ";";
                    outfile << std::fixed << std::setprecision(0) << flags_c[index] << ";";
                    outfile << std::fixed << std::setprecision(3) << score_c[index] << ";";
                    outfile << std::fixed << std::setprecision(0) << saved_c[index] << ";";
                    outfile << std::fixed << std::setprecision(0) << "-999;-999;-999;-999;-999;" << std::endl; 
                }
            }
        } // End CLUSTER PRESERVATION LOOP

        //---------------------------------------------------------------------
        // STRAY DATA REDEMPTION LOOP
        // Bringing back good observations that got caught in the wrong crowd
        int count_oi_r = 0;
        for(int curr=0; curr < s; curr++) {

            // break out if observation has NOT been flagged in the cluster preservation loop
            if(flags_c[curr] == 0) 
                continue;

            // get all neighbours that are close enough (consider only non-flagged neighbours)
            vec distances;
            ivec neighbour_indices = tree.get_neighbours_with_distance(lats[curr], lons[curr], outer_radius, distances, true);
            sct_with_fg_remove_flagged(neighbour_indices, distances, flags, flags_c, curr, num_max, 0, 0);
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
                int index = neighbour_indices[i];
                if(index == curr) { 
                    curr_in_box = i;
                    break;
                }
            }

            // call SCT with this box 
            vec lons_box = titanlib::subset(lons, neighbour_indices);
            vec elevs_box = titanlib::subset(elevs, neighbour_indices);
            vec lats_box = titanlib::subset(lats, neighbour_indices);
            vec values_box = titanlib::subset(values, neighbour_indices);
            vec backg_box = titanlib::subset(background_values, neighbour_indices);
            vec eps2_box = titanlib::subset(eps2, neighbour_indices);

            // perform statistical interpolation within the box
            vec ares(s_box,0);
            vec cvres(s_box,0);
            vec innov(s_box,0);
            float Dh_mean = 0.;
            sct_with_fg_oi_multi(ares, cvres, innov, Dh_mean, lons_box, lats_box, elevs_box, values_box, backg_box, eps2_box, vertical_scale, min_horizontal_scale, max_horizontal_scale, values_min, values_max, 1);
            count_oi_r++;
            // estimate observation error variance
            float sig2o = sct_with_fg_sig2_estimate(ares, innov, min_obs_var[curr], 0.5); 

            // save variables for diagnostics
            dh_r[curr] = Dh_mean;
            sig2_r[curr] = sig2o;

            // check the curr observation
            float score = cvres[curr_in_box] * ares[curr_in_box] / sig2o;
            assert(titanlib::is_valid(score));
            score_r[curr] = score;
            if((cvres[curr_in_box] < 0 && score > pos[curr]) || (cvres[curr_in_box] >= 0 && score > neg[curr])) {
                flags_r[curr] = 1;
            } else {
                saved_r[curr] = 1;
            }

            if(diagnostics) {
                for(int i = 0; i < s_box; i++) {
                    int index = neighbour_indices[i];
                    // "it;loop;curr;i;index;lon;lat;z;yo;yb;ya;yav;dh;sig2;flags_d;score_d;flags_c;score_c;saved_c;flags_r;score_r;saved_r;flags;sct_score;"
                    outfile << iteration << ";3;" << curr << ";" << i << ";" << index << ";";
                    outfile << std::fixed << std::setprecision(5) << lons_box[i] << ";" << lats_box[i] << ";";
                    outfile << std::fixed << std::setprecision(0) << elevs_box[i] << ";";
                    outfile << std::fixed << std::setprecision(2) << values_box[i] << ";" << backg_box[i] << ";";
                    outfile << std::fixed << std::setprecision(2) << ares[i] + values_box[i] << ";" << cvres[i] + values_box[i] << ";";
                    outfile << std::fixed << std::setprecision(0) << dh_r[curr] << ";";
                    outfile << std::fixed << std::setprecision(3) << sig2_r[curr] << ";";
                    outfile << std::fixed << std::setprecision(0) << "-999;";
                    outfile << std::fixed << std::setprecision(0) << flags_d[index] << ";";
                    outfile << std::fixed << std::setprecision(3) << score_d[index] << ";";
                    outfile << std::fixed << std::setprecision(0) << flags_c[index] << ";";
                    outfile << std::fixed << std::setprecision(3) << score_c[index] << ";";
                    outfile << std::fixed << std::setprecision(0) << saved_c[index] << ";";
                    outfile << std::fixed << std::setprecision(0) << flags_r[index] << ";";
                    outfile << std::fixed << std::setprecision(3) << score_r[index] << ";";
                    outfile << std::fixed << std::setprecision(0) << saved_r[index] << ";";
                    outfile << std::fixed << std::setprecision(0) << "-999;-999;" << std::endl; 
                }
            }

        } // End STRAY DATA REDEMPTION LOOP

        //---------------------------------------------------------------------
        // FLAG ASSIGNMENT STEP
        // Definitive flag is assigned to each observation based on the outcomes 
        // of the detection, cluster preservation, and redemption processes
        int tot_checked_d = 0;
        int tot_flagged_d = 0;
        int tot_flagged_c = 0;
        int tot_flagged_r = 0;
        int tot_saved_c = 0;
        int tot_saved_r = 0;
        int tot_thrown_out = 0;
        for(int curr=0; curr < s; curr++) {
            // Flagged is a previous SCT iteration
            if(flags[curr] != 0) {
                tot_thrown_out++;
                continue;
            }
            // Fresh new flagged observation
            if(flags_r[curr] == 1) {
                flags[curr] = 1;
                // Keep the highest score
                sct_score[curr] = std::max(sct_score[curr], score_r[curr]);
                thrown_out++;
                tot_thrown_out++;
            // Not flagged observation
            } else {
                // Keep the highest score
                sct_score[curr] = std::max(sct_score[curr], score_d[curr]);
            }
            // diagnostics
            if(diagnostics) {
                if(checked_d[curr]==1)
                    tot_checked_d++;
                if(flags_d[curr]==1)
                    tot_flagged_d++;
                if(flags_c[curr]==1)
                    tot_flagged_c++;
                if(flags_r[curr]==1)
                    tot_flagged_r++;
                if(saved_c[curr]==1)
                    tot_saved_c++;
                if(saved_r[curr]==1)
                    tot_saved_r++;
            }
        } // End FLAG ASSIGNMENT STEP

        if(diagnostics) {
            for(int curr=0; curr < s; curr++) {
                // "it;loop;curr;i;index;lon;lat;z;yo;yb;ya;yav;dh;sig2;flags_d;score_d;flags_c;score_c;saved_c;flags_r;score_r;saved_r;flags;sct_score;"
                outfile << iteration << ";4;" << curr << ";" << 0 << ";" << curr << ";";
                outfile << std::fixed << std::setprecision(5) << lons[curr] << ";" << lats[curr] << ";";
                outfile << std::fixed << std::setprecision(0) << elevs[curr] << ";";
                outfile << std::fixed << std::setprecision(2) << values[curr] << ";" << background_values[curr] << ";";
                outfile << std::fixed << std::setprecision(2) << "-999;-999;";
                outfile << std::fixed << std::setprecision(0) << "-999;";
                outfile << std::fixed << std::setprecision(3) << "-999;";
                outfile << std::fixed << std::setprecision(0) << "-999;";
                outfile << std::fixed << std::setprecision(0) << flags_d[curr] << ";";
                outfile << std::fixed << std::setprecision(3) << score_d[curr] << ";";
                outfile << std::fixed << std::setprecision(0) << flags_c[curr] << ";";
                outfile << std::fixed << std::setprecision(3) << score_c[curr] << ";";
                outfile << std::fixed << std::setprecision(0) << saved_c[curr] << ";";
                outfile << std::fixed << std::setprecision(0) << flags_r[curr] << ";";
                outfile << std::fixed << std::setprecision(3) << score_r[curr] << ";";
                outfile << std::fixed << std::setprecision(0) << saved_r[curr] << ";";
                outfile << std::fixed << std::setprecision(0) << flags[curr] << ";";
                outfile << std::fixed << std::setprecision(3) << sct_score[curr] << ";" << std::endl; 
            }
            std::cout << "==========================================" << std::endl;
            std::cout << "-- Detection Loop -- " << std::endl;
            std::cout << "# obs " << s << " Checked " << tot_checked_d << " Flagged " << tot_flagged_d << " Number of OI loop " << count_oi_d << std::endl;
            std::cout << "-- Cluster Preservation Loop -- " << std::endl;
            std::cout << "# checked " << tot_flagged_d << " Flagged " << tot_flagged_c << " Saved " << tot_saved_c << " Number of OI loop " << count_oi_c << std::endl;
            std::cout << "-- Stray Data Redemption Loop -- " << std::endl;
            std::cout << "# checked " << tot_flagged_c << " Flagged " << tot_flagged_r << " Saved " << tot_saved_r << " Number of OI loop " << count_oi_r << std::endl;
            std::cout << "-- Summary -- " << std::endl;
            std::cout << "Removing " << thrown_out << " (Total removed " << tot_thrown_out << " , " << std::setprecision(0) << 100*tot_thrown_out/s << "%)" << std::endl;
            double e_time0 = titanlib::clock();
            std::cout << "Time " << std::setprecision(0) << e_time0 - s_time0 << std::endl;
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

void sct_with_fg_remove_flagged(ivec& indices, vec& distances, const ivec& flags1, const ivec& flags2, const int& index_to_keep, const int& nmax, const int& f1, const int& f2) {

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

} // END sct_with_fg_remove_flagged

//
void sct_with_fg_oi(vec& ares, vec& cvres, vec& innov, float& Dh_mean, const vec& lons, const vec& lats, const vec& elevs, const vec& values, const vec& backg, const vec& eps2, const float& Dz, const float& Dh_min, const float& Dh_max, const float& vmin, const float& vmax) {

    int s_dim = values.size();

    boost::numeric::ublas::matrix<float> disth(s_dim, s_dim);
    boost::numeric::ublas::matrix<float> distz(s_dim, s_dim);
    boost::numeric::ublas::vector<float> Dh(s_dim);

    for(int i=0; i < s_dim; i++) {
        vec Dh_vector(s_dim);
        for(int j=0; j < s_dim; j++) {
            disth(i, j) = titanlib::calc_distance(lats[i], lons[i], lats[j], lons[j]);
            distz(i, j) = fabs(elevs[i] - elevs[j]);
            if(i != j) {
                if(i < j)
                    Dh_vector[j - 1] = disth(i, j);
                else if(i > j)
                    Dh_vector[j] = disth(i, j);
            }
        }
        Dh(i) = titanlib::compute_quantile(0.10, Dh_vector);
    }

    Dh_mean = std::accumulate(std::begin(Dh), std::end(Dh), 0.0) / Dh.size();
    if(Dh_mean < Dh_min) 
        Dh_mean = Dh_min;
    else if(Dh_mean > Dh_max)
        Dh_mean = Dh_max;

    boost::numeric::ublas::matrix<float> S(s_dim,s_dim);
    boost::numeric::ublas::matrix<float> Sinv(s_dim,s_dim);
    for(int i=0; i < s_dim; i++) {
        for(int j=0; j < s_dim; j++) {
            double corr_value = std::exp(-.5 * std::pow((disth(i, j) / Dh_mean), 2) - .5 * std::pow((distz(i, j) / Dz), 2));
            if(i==j) { // weight the diagonal?? (0.5 default)
                corr_value = corr_value + eps2[i];
            }
            S(i,j) = corr_value;
        }
    }

    boost::numeric::ublas::vector<float> d(s_dim);
    for(int i=0; i < s_dim; i++) {
        d(i) = values[i] - backg[i]; // innovation: difference between actual value and background
        innov[i] = values[i] - backg[i];
    }

    //   Beginning of real SCT
    bool b = titanlib::invert_matrix(S, Sinv);
//    if(b != true) {
//        // TODO: flag differently or give an error???
//        continue;
//    }

    // Unweight the diagonal
    for(int i=0; i < s_dim; i++) {
        S(i,i) -= eps2[i];
    }

    boost::numeric::ublas::vector<float> Zinv(s_dim), Sinv_d(s_dim), ainc(s_dim);

    Sinv_d = boost::numeric::ublas::prod(Sinv, d);

    ainc = boost::numeric::ublas::prod(S, Sinv_d); // analysis increment, ainc = ya - yb

    for(int i=0; i<s_dim; i++) {
        Zinv(i) = (1/Sinv(i,i)); 
        ares[i] = ainc(i)-d(i); // analysis residual, ares=ya-yb-(yo-yb) = ya - yo
        if((vmin!=vmax) && (ares[i] < (vmin-values[i])))
            ares[i] = vmin-values[i];
        if((vmin!=vmax) && (ares[i] > (vmax-values[i])))
            ares[i] = vmax-values[i];
    }

    for(int i=0; i<s_dim; i++) {
        cvres[i] = -1*Zinv(i) * Sinv_d(i);  // CVAres = yav - yo
        if((vmin!=vmax) && (cvres[i] < (vmin-values[i])))
            cvres[i] = vmin-values[i];
        if((vmin!=vmax) && (cvres[i] > (vmax-values[i])))
            cvres[i] = vmax-values[i];
    }

} // END sct_with_fg_oi

//
void sct_with_fg_oi_multi(vec& ares, vec& cvres, vec& innov, float& Dh_mean, const vec& lons, const vec& lats, const vec& elevs, const vec& values, const vec& backg, const vec& eps2, const float& Dz, const float& Dh_min, const float& Dh_max, const float& vmin, const float& vmax, const int& n_multi) {

    int s_dim = values.size();

    boost::numeric::ublas::matrix<float> disth(s_dim, s_dim);
    boost::numeric::ublas::matrix<float> distz(s_dim, s_dim);
    boost::numeric::ublas::vector<float> Dh(s_dim);

    for(int i=0; i < s_dim; i++) {
        vec Dh_vector(s_dim);
        for(int j=0; j < s_dim; j++) {
            disth(i, j) = titanlib::calc_distance(lats[i], lons[i], lats[j], lons[j]);
            distz(i, j) = fabs(elevs[i] - elevs[j]);
            if(i != j) {
                if(i < j)
                    Dh_vector[j - 1] = disth(i, j);
                else if(i > j)
                    Dh_vector[j] = disth(i, j);
            }
        }
        Dh(i) = titanlib::compute_quantile(0.10, Dh_vector);
    }

    Dh_mean = std::accumulate(std::begin(Dh), std::end(Dh), 0.0) / Dh.size();
    if(Dh_mean < Dh_min) 
        Dh_mean = Dh_min;
    else if(Dh_mean > Dh_max)
        Dh_mean = Dh_max;

    boost::numeric::ublas::matrix<float> S(s_dim,s_dim);
    boost::numeric::ublas::matrix<float> Sinv(s_dim,s_dim);
    for(int i=0; i < s_dim; i++) {
        for(int j=0; j < s_dim; j++) {
            double corr_value = std::exp(-.5 * std::pow((disth(i, j) / Dh_mean), 2) - .5 * std::pow((distz(i, j) / Dz), 2));
            if(i==j) { // weight the diagonal?? (0.5 default)
                corr_value = corr_value + eps2[i];
            }
            S(i,j) = corr_value;
        }
    }

    boost::numeric::ublas::vector<float> d(s_dim);
    for(int i=0; i < s_dim; i++) {
        d(i) = values[i] - backg[i]; // innovation: difference between actual value and background
        innov[i] = values[i] - backg[i];
    }

    //   Beginning of real SCT
    bool b = titanlib::invert_matrix(S, Sinv);
//    if(b != true) {
//        // TODO: flag differently or give an error???
//        continue;
//    }

    // Unweight the diagonal
    boost::numeric::ublas::vector<float> Zinv(s_dim);
    for(int i=0; i < s_dim; i++) {
        S(i,i) -= eps2[i];
        Zinv(i) = (1/Sinv(i,i)); 
    }

    // Multi loop
    for(int k=0; k < n_multi; k++) {

        boost::numeric::ublas::vector<float> Sinv_d(s_dim), ainc(s_dim);

        Sinv_d = boost::numeric::ublas::prod(Sinv, d);

        ainc = boost::numeric::ublas::prod(S, Sinv_d); // analysis increment, ainc = ya - yb

        for(int i=0; i < s_dim; i++) {
            // analysis residual, ares = ya-yb-(yo-yb) = ya - yo
            ares[i] = ainc(i)-d(i); 
            if((vmin!=vmax) && (ares[i] < (vmin-values[i])))
                ares[i] = vmin-values[i];
            if((vmin!=vmax) && (ares[i] > (vmax-values[i])))
                ares[i] = vmax-values[i];
            // cv-analysis residual, cvres = yav - yo
            cvres[i] = -1*Zinv(i) * Sinv_d(i);
            if((vmin!=vmax) && (cvres[i] < (vmin-values[i])))
                cvres[i] = vmin-values[i];
            if((vmin!=vmax) && (cvres[i] > (vmax-values[i])))
                cvres[i] = vmax-values[i];
            // update innovation: difference between actual value and cv-value
            d(i) = -1 * cvres[i]; 
        }
    }  // End multi loop

} // END sct_with_fg_oi

//
float sct_with_fg_sig2_estimate(const vec& ares, const vec& innov, const float& vmin, const float& q) {

    int s_dim = ares.size();

    std::vector<float> sig2_temp(s_dim);
    for(int i=0; i<s_dim; i++) 
        sig2_temp[i] = -1*ares[i]*innov[i];
 
    float sig2 = titanlib::compute_quantile( q, sig2_temp);
    if(sig2 < vmin) 
        sig2 = vmin;
    return sig2;
} // END sct_with_fg_sig2_estimate

//
float sct_with_fg_sig2_estimate_alt(const vec& values, const float& vmin) {

    // pseudo-variance
//    float sig2 = std::pow( ((titanlib::compute_quantile( 0.75, values) - titanlib::compute_quantile( 0.25, values)) / 1.349), 2);
    float sig2 = std::pow( titanlib::compute_quantile( 0.9, values), 2);

    if(sig2 <vmin) 
        sig2 = vmin;

    return sig2;
} // END sct_with_fg_sig2_estimate


