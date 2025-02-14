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
// for the debug file
#include <iostream>
#include <fstream>  // For file handling
//#include <filesystem>  // C++17 for checking file existence and deleting
#include <cstdio>      // For std::remove

using namespace titanlib;

// helpers
//void remove_flagged_copy(ivec& indices, vec& distances, const ivec& flags);
void remove_flagged_copy(ivec& indices, vec& distances, const ivec& flags, const int index_to_keep);

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
        float vertical_scale,
        const vec& pos,
        const vec& neg,
        const vec& eps2,
        const vec& min_obs_var,
        vec& gross_error_score,
        const ivec& obs_to_check) {

    //debug
    bool debug = true;
    std::ofstream outfile; 
    if (debug) {
        const std::string filename = "/home/cristianl/test.txt";
        // Check if the file exists and delete it
        if (std::ifstream(filename)) {
            std::remove(filename.c_str());
            std::cout << "File existed and was deleted.\n";
        }
        // Check if the file exists and delete it
//        if (std::filesystem::exists(filename)) {
//            std::filesystem::remove(filename);
//            std::cout << "File existed and was deleted.\n";
//        }
        outfile.open(filename, std::ios_base::app); // Open file in append mode
        if (!outfile) { // Check if the file opened successfully
            std::cerr << "Error opening file!" << std::endl;
        }
        outfile << "it;loop;curr;i;index;lon;lat;z;yo;yb;ares;cvres;sig2_1st;flags_1st;ges_1st;dh;sig2_2nd;flags_2nd;ges_2nd;t_ges2nd;saved;" << std::endl;
    }

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
    gross_error_score.clear();
    gross_error_score.resize(s, 0);

    titanlib::KDTree tree(lats, lons);

    // Flag stations without elevation
    for(int curr=0; curr < s; curr++) {
        if(!titanlib::is_valid(elevs[curr])) {
            flags[curr] = 1;
        }
    }
    for(int iteration = 0; iteration < num_iterations; iteration++) {
        double s_time0 = titanlib::clock();

        int thrown_out = 0; // reset this number each loop (this is for breaking if we don't throw anything new out)
        int saved = 0; // reset this number each loop (this is for knowing how many observations have been given a 2nd chance)

        ivec checked(s, 0);  // Keep track of which observations have been checked
        ivec flags_1st(s, 0);  // Keep track of which observations have been checked in this iteration only
        vec ges_1st(s, 0); // GES in the screening loop
        vec ges_2nd(s, 0); // GES in the second chance loop
        // debug
        vec dh_deb(s, 0); 
        vec sig2_1st_deb(s, 0); 
        vec sig2_2nd_deb(s, 0); 
        vec saved_deb(s, 0); 

        //---------------------------------------------------------------------
        // SCREENING LOOP
        int count_oi_1st = 0;
        for(int curr=0; curr < s; curr++) {
            if(obs_to_check.size() == s && obs_to_check[curr] != 1) {
                checked[curr] = 1;
                continue;
            }

            // break out if station already flagged
            if(flags[curr] != 0) {
                checked[curr] = 1;
                continue;
            }
            if(checked[curr] > 0) {
                continue;
            }
            
            // break out if observation close to first guess
            float buffer = 2 * std::sqrt(min_obs_var[curr]);
            if((values[curr] >= (background_values[curr]-buffer)) && (values[curr] <= (background_values[curr]+buffer))) {
                checked[curr] = 1;
                continue;
            }

            // get all neighbours that are close enough
            vec distances;
            ivec neighbour_indices = tree.get_neighbours_with_distance(lats[curr], lons[curr], outer_radius, distances, true);
            remove_flagged_copy(neighbour_indices, distances, flags, -1);

            if(neighbour_indices.size() > num_max) {

                int N = neighbour_indices.size();
                std::vector<std::pair<float,int> > pairs(N);
                for(int i = 0; i < neighbour_indices.size(); i++) {
                    pairs[i] = std::pair<float, int>(distances[i], neighbour_indices[i]);
                }
                std::sort(pairs.begin(), pairs.end(), titanlib::sort_pair_first<float,int>());
                distances.clear();
                neighbour_indices.clear();
                distances.resize(num_max);
                neighbour_indices.resize(num_max);
                for(int i = 0; i < num_max; i++) {
                    distances[i] = pairs[i].first;
                    neighbour_indices[i] = pairs[i].second;
                }
            }
            int s_box = neighbour_indices.size();
            if(s_box < num_min) {
                checked[curr] = 1;
                // flag as isolated? 
                continue; // go to next station, skip this one
            }


            // call SCT with this box 
            vec lons_box = titanlib::subset(lons, neighbour_indices);
            vec elevs_box = titanlib::subset(elevs, neighbour_indices);
            vec lats_box = titanlib::subset(lats, neighbour_indices);
            vec values_box = titanlib::subset(values, neighbour_indices);
            vec backg_box = titanlib::subset(background_values, neighbour_indices);
            vec eps2_box = titanlib::subset(eps2, neighbour_indices);
            // the thing to flag is at "curr", ano not included in the box

            boost::numeric::ublas::matrix<float> disth(s_box, s_box);
            boost::numeric::ublas::matrix<float> distz(s_box, s_box);
            boost::numeric::ublas::vector<float> Dh(s_box);

            for(int i=0; i < s_box; i++) {
                vec Dh_vector(s_box);
                for(int j=0; j < s_box; j++) {
                    disth(i, j) = titanlib::calc_distance(lats_box[i], lons_box[i], lats_box[j], lons_box[j]);
                    distz(i, j) = fabs(elevs_box[i] - elevs_box[j]);
                    if(i != j) {
                        if(i < j)
                            Dh_vector[j - 1] = disth(i, j);
                        else if(i > j)
                            Dh_vector[j] = disth(i, j);
                    }
                }
                Dh(i) = titanlib::compute_quantile(0.10, Dh_vector);
            }

            double Dh_mean = std::accumulate(std::begin(Dh), std::end(Dh), 0.0) / Dh.size();
            dh_deb[curr] = Dh_mean;
            if(Dh_mean < min_horizontal_scale) {
                Dh_mean = min_horizontal_scale;
            }

            boost::numeric::ublas::matrix<float> S(s_box,s_box);
            boost::numeric::ublas::matrix<float> Sinv(s_box,s_box);
            for(int i=0; i < s_box; i++) {
                for(int j=0; j < s_box; j++) {
                    double value = std::exp(-.5 * std::pow((disth(i, j) / Dh_mean), 2) - .5 * std::pow((distz(i, j) / vertical_scale), 2));
                    if(i==j) { // weight the diagonal?? (0.5 default)
                        value = value + eps2_box[i];
                    }
                    S(i,j) = value;
                }
            }

            boost::numeric::ublas::vector<float> d(s_box);
            for(int i=0; i < s_box; i++) {
                d(i) = values_box[i] - backg_box[i]; // innovation: difference between actual value and background
            }

            //   Beginning of real SCT
            bool b = titanlib::invert_matrix(S, Sinv);
            if(b != true) {
                // TODO: flag differently or give an error???
                continue;
            }

            // Unweight the diagonal
            for(int i=0; i < s_box; i++) {
                S(i,i) -= eps2_box[i];
            }

            boost::numeric::ublas::vector<float> Zinv(s_box), Sinv_d(s_box), ainc(s_box), ares(s_box);

            Sinv_d = boost::numeric::ublas::prod(Sinv, d);

            ainc = boost::numeric::ublas::prod(S, Sinv_d); // analysis increment, ainc = ya - yb

            for(int i=0; i<s_box; i++) {
                Zinv(i) = (1/Sinv(i,i)); 
                ares(i) = ainc(i)-d(i); // analysis residual, ares=ya-yb-(yo-yb) = ya - yo
                if((values_min!=values_max) && (ares(i) < (values_min-values_box[i])))
                    ares(i) = values_min-values_box[i];
                if((values_min!=values_max) && (ares(i) > (values_max-values_box[i])))
                    ares(i) = values_max-values_box[i];
            }

            boost::numeric::ublas::vector<float> cvres(s_box);
            for(int i=0; i<s_box; i++) {
                cvres(i) = -1*Zinv(i) * Sinv_d(i);  // CVAres = yav - yo
                if((values_min!=values_max) && (cvres(i) < (values_min-values_box[i])))
                    cvres(i) = values_min-values_box[i];
                if((values_min!=values_max) && (cvres(i) > (values_max-values_box[i])))
                    cvres(i) = values_max-values_box[i];
            }

            std::vector<float> sig2o_temp(s_box);
            for(int i=0; i<s_box; i++) {
                sig2o_temp[i] = -1*ares(i)*d(i);
            }
 
            float sig2o = titanlib::compute_quantile( 0.5, sig2o_temp);
            sig2_1st_deb[curr] = sig2o;
            if(sig2o < min_obs_var[curr]) {
                sig2o = min_obs_var[curr];
            }

            int ccount = 0;
            for(int i = 0; i < s_box; i++) {
                int index = neighbour_indices[i];
                if(obs_to_check.size() == s && obs_to_check[index] != 1) {
                    checked[curr] = 1;
                    continue;
                }
                float dist = distances[i];
                if(dist <= inner_radius) {
                    float ges = cvres(i) * ares(i) / sig2o;
                    assert(titanlib::is_valid(ges));
                    gross_error_score[index] = std::max(ges, gross_error_score[index]);
                    ges_1st[index] = ges;
                    // condition to identify a candidate for flagging
                    if((cvres(i) < 0 && ges > pos[index]) || (cvres(i) >= 0 && ges > neg[index])) {
                        flags_1st[index] = 1;
                    }
                    checked[index] = 1;
                    ccount++;
                }
            }
            count_oi_1st++;

            // debugging
            if(debug) {
                for(int i = 0; i < s_box; i++) {
                    int index = neighbour_indices[i];
                    outfile << iteration << ";1;" << curr << ";" << i << ";" << index << ";";
                    outfile << std::fixed << std::setprecision(5) << lons_box[i] << ";" << lats_box[i] << ";";
                    outfile << std::fixed << std::setprecision(0) << elevs_box[i] << ";";
                    outfile << std::fixed << std::setprecision(2) << values_box[i] << ";" << backg_box[i] << ";";
                    outfile << std::fixed << std::setprecision(2) << ares(i) << ";" << cvres(i) << ";";
                    outfile << std::fixed << std::setprecision(3) << sig2_1st_deb[curr] << ";";
                    outfile << std::fixed << std::setprecision(0) << flags_1st[index] << ";";
                    outfile << std::fixed << std::setprecision(3) << ges_1st[index] << ";";
                    outfile << std::fixed << std::setprecision(0) << dh_deb[curr] << ";0;0;0;0;0;" << std::endl; 
                }
            }

        } // End screening loop over observations to check

        //---------------------------------------------------------------------
        // provisional flagging of all suspect observations, such that they are excluded in the next remove_flagged_copy
        for(int curr=0; curr < s; curr++) {
          if (flags_1st[curr]==1) 
              flags[curr] = 1;
        }

        //---------------------------------------------------------------------
        // SECOND CHANCE
        // save suspect observations if their GESs is not too different from the neighbours
        int count_oi_2nd = 0;
        for(int curr=0; curr < s; curr++) {
            if(flags_1st[curr] != 1) {
                continue;
            }

            // get all neighbours that are close enough
            vec distances;
            ivec neighbour_indices = tree.get_neighbours_with_distance(lats[curr], lons[curr], outer_radius, distances, true);
            remove_flagged_copy(neighbour_indices, distances, flags, curr);

            if(neighbour_indices.size() > num_max) {

                int N = neighbour_indices.size();
                std::vector<std::pair<float,int> > pairs(N);
                for(int i = 0; i < neighbour_indices.size(); i++) {
                    pairs[i] = std::pair<float, int>(distances[i], neighbour_indices[i]);
                }
                std::sort(pairs.begin(), pairs.end(), titanlib::sort_pair_first<float,int>());
                distances.clear();
                neighbour_indices.clear();
                distances.resize(num_max);
                neighbour_indices.resize(num_max);
                for(int i = 0; i < num_max; i++) {
                    distances[i] = pairs[i].first;
                    neighbour_indices[i] = pairs[i].second;
                }
            }

            // call SCT with this box 
            vec lons_box = titanlib::subset(lons, neighbour_indices);
            vec elevs_box = titanlib::subset(elevs, neighbour_indices);
            vec lats_box = titanlib::subset(lats, neighbour_indices);
            vec ges_box = titanlib::subset(ges_1st, neighbour_indices);
            vec eps2_box = titanlib::subset(eps2, neighbour_indices);
            // just used for debugging (next 2 lines)
            vec values_box = titanlib::subset(values, neighbour_indices);
            vec backg_box = titanlib::subset(background_values, neighbour_indices);
            int s_box = neighbour_indices.size();
            // the thing to flag is at "curr", ano not included in the box

            boost::numeric::ublas::matrix<float> disth(s_box, s_box);
            boost::numeric::ublas::matrix<float> distz(s_box, s_box);
            boost::numeric::ublas::vector<float> Dh(s_box);

            for(int i=0; i < s_box; i++) {
                vec Dh_vector(s_box);
                for(int j=0; j < s_box; j++) {
                    disth(i, j) = titanlib::calc_distance(lats_box[i], lons_box[i], lats_box[j], lons_box[j]);
                    distz(i, j) = fabs(elevs_box[i] - elevs_box[j]);
                    if(i != j) {
                        if(i < j)
                            Dh_vector[j - 1] = disth(i, j);
                        else if(i > j)
                            Dh_vector[j] = disth(i, j);
                    }
                }
                Dh(i) = titanlib::compute_quantile(0.10, Dh_vector);
            }

            double Dh_mean = std::accumulate(std::begin(Dh), std::end(Dh), 0.0) / Dh.size();
            if(Dh_mean < min_horizontal_scale) {
                Dh_mean = min_horizontal_scale;
            }

            boost::numeric::ublas::matrix<float> S(s_box,s_box);
            boost::numeric::ublas::matrix<float> Sinv(s_box,s_box);
            for(int i=0; i < s_box; i++) {
                for(int j=0; j < s_box; j++) {
                    double value = std::exp(-.5 * std::pow((disth(i, j) / Dh_mean), 2) - .5 * std::pow((distz(i, j) / vertical_scale), 2));
                    if(i==j) { // weight the diagonal?? (0.5 default)
                        value = value + eps2_box[i];
                    }
                    S(i,j) = value;
                }
            }

            boost::numeric::ublas::vector<float> ges_tmp(s_box);
            for(int i=0; i < s_box; i++) {
                ges_tmp(i) = ges_box[i]; // innovation: difference between actual value and background
            }
            
            //   Beginning of real SCT
            bool b = titanlib::invert_matrix(S, Sinv);
            if(b != true) {
                // TODO: flag differently or give an error???
                continue;
            }

            // Unweight the diagonal
            for(int i=0; i < s_box; i++) {
                S(i,i) -= eps2_box[i];
            }

            boost::numeric::ublas::vector<float> Zinv(s_box), Sinv_ges(s_box), analysis(s_box), ares(s_box);

            Sinv_ges = boost::numeric::ublas::prod(Sinv, ges_tmp);

            analysis = boost::numeric::ublas::prod(S, Sinv_ges); 

            for(int i=0; i<s_box; i++) {
                Zinv(i) = (1/Sinv(i,i)); 
                ares(i) = analysis(i)-ges_box[i]; 
            }

            boost::numeric::ublas::vector<float> cvres(s_box);
            for(int i=0; i<s_box; i++) {
                cvres(i) = -1*Zinv(i) * Sinv_ges(i);
            }

            // pseudo-variance
            double sig2_2nd = std::pow( ((titanlib::compute_quantile( 0.75, ges_box) - titanlib::compute_quantile( 0.25, ges_box)) / 1.349), 2);
            sig2_2nd_deb[curr] = sig2_2nd;

            if(sig2_2nd < 0.1) {
                sig2_2nd = 0.1;
            }

            int ccount = 0;
            double alpha = 0.05; // Significance level (5%)
            // remember that "num_min must be > 1" and s_box is greather or euqal to num_min
            double df = s_box-1;
            // Create a Student's t-distribution
            boost::math::students_t tdistro(df);
            // Compute critical t-value for two-tailed test (alpha/2 for both tails)
            double t_critical = boost::math::quantile(boost::math::complement(tdistro, alpha / 2));
            for(int i = 0; i < s_box; i++) {
                int index = neighbour_indices[i];
                if(index == curr) {
                    float dist = distances[i];
                    // Compute test statistic (t-score)
                    float ges = cvres(i) * ares(i) / sig2_2nd;
                    assert(titanlib::is_valid(ges));
                    ges_2nd[index] = ges;
                    // float p_value = 2 * boost::math::cdf(boost::math::complement(tdistro, std::fabs(ges)));
                    // Reject the null hypothesis: Gross Error score is significantly different from its expected value
                    if(std::fabs(ges) > t_critical) {
//                        flags[index] = 1;
                        thrown_out++;
                    // Fail to reject the null hypothesis: No significant difference in Gross Error score
                    } else {
                        flags[index] = 0;
                        saved_deb[index] = 1;
                        saved++;
                    }
                }
            }
            count_oi_2nd++;

            // debugging
//            for(int i = 0; i < s_box; i++) {
//              int index = neighbour_indices[i];
//              std::cout << iteration << ";2;" << curr << ";" << i << ";" << lons_box[i] << ";" << lats_box[i] << ";" << elevs_box[i] << ";" << values_box[i] << ";" << backg_box[i] << ";" << ares(i) << ";" << cvres(i) << ";0;0;0;0;" << sig2_2nd_deb[curr] << ";" << flags[index] << ";" << ges_2nd[index] << ";" << t_critical << ";" << std::endl;
//            }
            // debugging
            if(debug) {
                for(int i = 0; i < s_box; i++) {
                    int index = neighbour_indices[i];
                    outfile << iteration << ";2;" << curr << ";" << i << ";" << index << ";";
                    outfile << std::fixed << std::setprecision(5) << lons_box[i] << ";" << lats_box[i] << ";";
                    outfile << std::fixed << std::setprecision(0) << elevs_box[i] << ";";
                    outfile << std::fixed << std::setprecision(2) << values_box[i] << ";" << backg_box[i] << ";";
                    outfile << std::fixed << std::setprecision(2) << ares(i) << ";" << cvres(i) << ";0;0;0;0;";
                    outfile << std::fixed << std::setprecision(3) << sig2_2nd_deb[curr] << ";";
                    outfile << std::fixed << std::setprecision(0) << flags[index] << ";";
                    outfile << std::fixed << std::setprecision(3) << ges_2nd[index] << ";" << t_critical << ";";
                    outfile << std::fixed << std::setprecision(0) << saved_deb[index] << ";" << std::endl;
                }
            }

        } // End second chance loop over observations to check 
        if(thrown_out == 0) {
            if(iteration + 1 < num_iterations) {
//                std::cout << "Stopping early after " << iteration + 1<< " iterations" << std::endl;
            }
            break;
        }
        std::cout << "Removing " << thrown_out << " Saving " << saved << " Number of OI 1st loop " << count_oi_1st << " Number of OI 2nd loop " << count_oi_2nd << std::endl;
//        double e_time0 = titanlib::clock();
//        std::cout << "Time " << e_time0 - s_time0 << std::endl;
    } // End loop over iterations

    return flags;
}
// end SCT //

//----------------------------------------------------------------------------//
// HELPER FUNCTIONS

void remove_flagged_copy(ivec& indices, vec& distances, const ivec& flags, const int index_to_keep) {
    ivec indices_new;
    vec distances_new;
    indices_new.reserve(indices.size());
    distances_new.reserve(indices.size());
    for(int i=0; i<indices.size(); i++) {
        if(flags[indices[i]] == 0 || indices[i] == index_to_keep) {
            indices_new.push_back(indices[i]);
            distances_new.push_back(distances[i]);
        }
    }
    indices = indices_new;
    distances = distances_new;
}
