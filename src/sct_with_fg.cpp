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

using namespace titanlib;

// helpers
void remove_flagged_copy(ivec& indices, vec& distances, const ivec& flags);

// start SCT //
ivec titanlib::sct_with_fg(const Points& points,
        const vec& values,
        const vec& background_values,
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
        vec& rep, // Can we remove this?
        const ivec& obs_to_check) {

    //debug
//    std::ofstream outfile;
//    outfile.open("/home/cristianl/test.txt", std::ios_base::app);
//    outfile << "Data";

    std::cout << "number of observations " << std::endl;
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
    rep.clear();
    rep.resize(s, 0);

    titanlib::KDTree tree(lats, lons);

    // Flag stations without elevation
    for(int curr=0; curr < s; curr++) {
        if(!titanlib::is_valid(elevs[curr])) {
            flags[curr] = 1;
        }
    }
    for(int iteration = 0; iteration < num_iterations; iteration++) {
        std::cout << "================= ITERATION ===== " << iteration << " ========================" << std::endl;
        std::cout << "number of observations " << s << std::endl;
        double s_time0 = titanlib::clock();

        int thrown_out = 0; // reset this number each loop (this is for breaking if we don't throw anything new out)
        int saved = 0; // reset this number each loop (this is for knowing how many observations have been given a 2nd chance)

        ivec checked(s, 0);  // Keep track of which observations have been checked
        ivec flags_tmp(s, 0);  // Keep track of which observations have been checked in this iteration only
        vec ges_tmp(s); // GES in this iteration only

        int count_oi = 0;
        for(int curr=0; curr < s; curr++) {
            std::cout << "----------------- OBSERVATION ----- " << curr << " ----------------------" << std::endl;
            std::cout << "observation background = " << values[curr] << " " << background_values[curr] << std::endl;
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
            if((values[curr] >= (background_values[curr]-min_obs_var[curr])) && (values[curr] <= (background_values[curr]+min_obs_var[curr]))) {
                checked[curr] = 1;
                continue;
            }

            // get all neighbours that are close enough
            vec distances;
            ivec neighbour_indices = tree.get_neighbours_with_distance(lats[curr], lons[curr], outer_radius, distances, true);
            remove_flagged_copy(neighbour_indices, distances, flags);

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
            std::cout << "s_box = " << s_box << std::endl;
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
            std::cout << "Dh_mean = " << Dh_mean << std::endl;
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

            boost::numeric::ublas::vector<float> Zinv(s_box), Sinv_d(s_box), ares_temp(s_box), ares(s_box);

            Sinv_d = boost::numeric::ublas::prod(Sinv, d);

            ares_temp = boost::numeric::ublas::prod(S, Sinv_d); // ares = ya - yb

            for(int i=0; i<s_box; i++) {
                Zinv(i) = (1/Sinv(i,i)); 
                ares(i) = ares_temp(i)-d(i); // Ares=ya-yb-(yo-yb) = ya - yo
            }

            boost::numeric::ublas::vector<float> cvres(s_box);
            vec cvanalysis(s_box);
            for(int i=0; i<s_box; i++) {
                cvres(i) = -1*Zinv(i) * Sinv_d(i);  // CVAres = yav - yo
                cvanalysis[i] = cvres(i) + values_box[i];
                std::cout << "outer circle, i observation background analysis cvanalysis = " << i << " " << values_box[i] << " " << backg_box[i] << " " << ares(i) + values_box[i] << " " << cvanalysis[i] << std::endl;
            }
            double sig2obs = std::pow( ((titanlib::compute_quantile( 0.75, values_box) - titanlib::compute_quantile( 0.25, values_box)) / 1.349), 2);
            double sig2cva = std::pow( ((titanlib::compute_quantile( 0.75, cvanalysis) - titanlib::compute_quantile( 0.25, cvanalysis)) / 1.349), 2);
            std::cout << "sig2obs sig2cva sig2sum = " << sig2obs << " " << sig2cva << " " << sig2obs+sig2cva << std::endl;

            double sig2o = 0;
            boost::numeric::ublas::vector<float> sig2o_temp(s_box), negAres_temp(s_box);
            for(int i=0; i<s_box; i++) {
                negAres_temp(i)=-1*ares(i);
                sig2o_temp(i) = d(i)*negAres_temp(i);
                sig2o += sig2o_temp(i);
            }
 
            // Convert to std::vector<float>
            std::vector<float> sig2o_temp1(sig2o_temp.begin(), sig2o_temp.end());

//            sig2o = sig2o/s_box;
//            std::cout << "sig2o = " << sig2o << std::endl;
            sig2o = titanlib::compute_quantile( 0.5, sig2o_temp1);
            std::cout << "sig2o = " << sig2o << std::endl;
            if(sig2o < min_obs_var[curr]) {
                sig2o = min_obs_var[curr];
            }

            // boost::numeric::ublas::vector<float> ges(s_box);
            // for(int i=0; i<s_box; i++) {
            //     ges(i) = cvres(i)*ares(i) / sig2o;
            // }
            int ccount = 0;
            for(int i = 0; i < s_box; i++) {
                int index = neighbour_indices[i];
                if(obs_to_check.size() == s && obs_to_check[index] != 1) {
                    checked[curr] = 1;
                    continue;
                }
                float dist = distances[i];
                if(dist <= inner_radius) {
                    std::cout << "................... inner circle ..... " << i << " " << index << " ......................." << std::endl;
                    float ges = cvres(i) * ares(i) / sig2o;
                    std::cout << "ges= " << ges << std::endl;
                    assert(titanlib::is_valid(ges));
                    gross_error_score[index] = std::max(ges, gross_error_score[index]);
                    ges_tmp[index] = ges;
                    if((cvres(i) < 0 && ges > pos[index]) || (cvres(i) >= 0 && ges > neg[index])) {
//                        flags[index] = 1;
                        flags_tmp[index] = 1;
//                        thrown_out++;
                    }
                    std::cout << "observation background ges flags_tmp= " << values_box[i] << " " << backg_box[i] << " " << ges << " " << flags_tmp[index] << std::endl;
                    checked[index] = 1;
                    ccount++;
                }
            }
            count_oi++;
        } // End loop over observations to check

        // SECOND CHANCE
        // save suspect observations if their GESs is not too different from the neighbours
        for(int curr=0; curr < s; curr++) {
            if(flags_tmp[curr] != 1) {
                continue;
            }
            std::cout << "----------------- 2ND CHANCE OBSERVATION ----- " << curr << " ----------------------" << std::endl;
            std::cout << "observation background = " << values[curr] << " " << background_values[curr] << std::endl;

            // get all neighbours that are close enough
            vec distances;
            ivec neighbour_indices = tree.get_neighbours_with_distance(lats[curr], lons[curr], outer_radius, distances, true);
            remove_flagged_copy(neighbour_indices, distances, flags);

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
            vec ges_box = titanlib::subset(ges_tmp, neighbour_indices);
            vec eps2_box = titanlib::subset(eps2, neighbour_indices);
            // just used for debugging (next 2 lines)
            vec values_box = titanlib::subset(values, neighbour_indices);
            vec backg_box = titanlib::subset(background_values, neighbour_indices);
            int s_box = neighbour_indices.size();
            std::cout << "s_box = " << s_box << std::endl;
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
            std::cout << "Dh_mean = " << Dh_mean << std::endl;
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
            
            // debug
            boost::numeric::ublas::vector<float> ges_box1(s_box);
            for(int i=0; i < s_box; i++) {
                std::cout << "outer circle, i obs backg ges = " << i << " " << values_box[i] << " " << backg_box[i] << " " << ges_box[i] << std::endl;
                ges_box1(i) = ges_box[i];
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

            Sinv_ges = boost::numeric::ublas::prod(Sinv, ges);

            analysis = boost::numeric::ublas::prod(S, Sinv_ges); 

//            for(int i=0; i<s_box; i++) {
//                double acc = 0;
//                for(int j=0; j<s_box; j++) {
//                    acc += Sinv(i,j)*ges_box1(j);
//                }
//                Sinv_ges(i) = acc;
//            }
//            for(int i=0; i<s_box; i++) {
//                double acc = 0;
//                for(int j=0; j<s_box; j++) {
//                    acc += S(i,j)*Sinv_ges(j);
//                }
//                analysis(i) = acc;
//            }
            for(int i=0; i<s_box; i++) {
                Zinv(i) = (1/Sinv(i,i)); 
                ares(i) = analysis(i)-ges_box1(i); 
            }

            boost::numeric::ublas::vector<float> cvres(s_box);
            for(int i=0; i<s_box; i++) {
                cvres(i) = -1*Zinv(i) * Sinv_ges(i);
            }

            // Convert to std::vector<float>
            double sig2_2nd = std::pow( ((titanlib::compute_quantile( 0.75, ges_box) - titanlib::compute_quantile( 0.25, ges_box)) / 1.349), 2);
            std::cout << "sig2_2nd = " << sig2_2nd << std::endl;

            if(sig2_2nd < 0.1) {
                sig2_2nd = 0.1;
            }

            // boost::numeric::ublas::vector<float> ges(s_box);
            // for(int i=0; i<s_box; i++) {
            //     ges(i) = cvres(i)*ares(i) / sig2_2nd;
            // }

            int ccount = 0;
            double alpha = 0.05; // Significance level (5%)
            for(int i = 0; i < s_box; i++) {
                int index = neighbour_indices[i];
                float dist = distances[i];
                double df = s_box-1;
                // Create a Student's t-distribution
                boost::math::students_t tdistro(df);
                // Compute critical t-value for two-tailed test (alpha/2 for both tails)
                double t_critical = boost::math::quantile(boost::math::complement(tdistro, alpha / 2));
                if((dist <= inner_radius) && (flags_tmp[index] == 1)) {
                    std::cout << "................... inner circle ..... " << i << " " << index << " ......................." << std::endl;
                    // Compute test statistic (t-score)
                    float ges_2nd = cvres(i) * ares(i) / sig2_2nd;
                    assert(titanlib::is_valid(ges_2nd));
//                    float p_value = 2 * boost::math::cdf(boost::math::complement(tdistro, std::fabs(ges_2nd)));
                    // Reject the null hypothesis: Gross Error score is significantly different
                    if(std::fabs(ges_2nd) > t_critical) {
                        flags[index] = 1;
                        thrown_out++;
                    // Fail to reject the null hypothesis: No significant difference in Gross Error score
                    } else {
                        saved++;
                        std::cout << "***" << std::endl;
                    }
                    std::cout << "cvres ares= " << cvres[i] << " " << ares[i] << std::endl;
                    std::cout << "observation background ges_2nd flags= " << values[index] << " " << background_values[index] << " " << ges_2nd << " " << flags[index] << std::endl;
                }
            }
        } // End loop over observations to check 
        if(thrown_out == 0) {
            if(iteration + 1 < num_iterations) {
                std::cout << "Stopping early after " << iteration + 1<< " iterations" << std::endl;
            }
            break;
        }
        std::cout << "Removing " << thrown_out << " Saving " << saved << " Number of OI " << count_oi << std::endl;
        double e_time0 = titanlib::clock();
        // std::cout << e_time0 - s_time0 << std::endl;
    } // End loop over iterations

    return flags;
}
// end SCT //

//----------------------------------------------------------------------------//
// HELPER FUNCTIONS

void remove_flagged_copy(ivec& indices, vec& distances, const ivec& flags) {
    ivec indices_new;
    vec distances_new;
    indices_new.reserve(indices.size());
    distances_new.reserve(indices.size());
    for(int i=0; i<indices.size(); i++) {
        if(flags[indices[i]] == 0 ) {
            indices_new.push_back(indices[i]);
            distances_new.push_back(distances[i]);
        }
    }
    indices = indices_new;
    distances = distances_new;
}
