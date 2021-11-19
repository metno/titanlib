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

using namespace titanlib;

// helpers
void remove_flagged(ivec& indices, vec& distances, const ivec& flags);

// start SCT //
ivec titanlib::sct(const Points& points,
        const vec& values,
        // determine if we have too many or too few observations
        // (too many means we can reduce the distance, too few mean isolation problem and cannot flag?)
        int num_min,
        int num_max,
        // first find everything close to the point that we are testing (maxdist)
        float inner_radius,
        float outer_radius,
        int num_iterations,
        int num_min_prof,
        float min_elev_diff,
        float min_horizontal_scale,
        float vertical_scale,
        const vec& pos,
        const vec& neg,
        const vec& eps2,
        vec& prob_gross_error,
        vec& rep,
        const ivec& obs_to_check) {

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
    if(obs_to_check.size() > 0 && obs_to_check.size() != s)
        throw std::invalid_argument("obs_to_check must empty or have the same size as values");
    if(num_min < 2)
        throw std::invalid_argument("num_min must be > 1");
    if(num_max < num_min)
        throw std::invalid_argument("num_max must be > num_min");
    if(num_iterations < 1)
        throw std::invalid_argument("num_iterations must be >= 1");
    if(min_elev_diff <= 0)
        throw std::invalid_argument("min_elev_diff must be > 0");
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
        if(pos[i] < 0)
            throw std::invalid_argument("All pos values must be >= 0");
        if(neg[i] < 0)
            throw std::invalid_argument("All neg values must be >= 0");
    }

    ivec flags(s, 0);
    prob_gross_error.clear();
    prob_gross_error.resize(s, 0);
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
        double s_time0 = titanlib::clock();

        int thrown_out = 0; // reset this number each loop (this is for breaking if we don't throw anything new out)

        ivec checked(s, 0);  // Keep track of which observations have been checked
        int count_oi = 0;
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

            // get all neighbours that are close enough
            vec distances;
            ivec neighbour_indices = tree.get_neighbours_with_distance(lats[curr], lons[curr], outer_radius, distances, true);
            remove_flagged(neighbour_indices, distances, flags);

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
            if(neighbour_indices.size() < num_min) {
                checked[curr] = 1;
                // flag as isolated? 
                continue; // go to next station, skip this one
            }

            // call SCT with this box 
            vec lons_box = titanlib::subset(lons, neighbour_indices);
            vec elevs_box = titanlib::subset(elevs, neighbour_indices);
            vec lats_box = titanlib::subset(lats, neighbour_indices);
            vec values_box = titanlib::subset(values, neighbour_indices);
            vec eps2_box = titanlib::subset(eps2, neighbour_indices);
            int s_box = neighbour_indices.size();
            // the thing to flag is at "curr", ano not included in the box

            // Compute the background
            vec vp;
            if(num_min_prof >= 0) {
                vp = titanlib::background(elevs_box, values_box, num_min_prof, min_elev_diff, titanlib::MV, titanlib::MV, titanlib::VerticalProfile, vec(), ivec(), false);
            }
            else {
                double meanT = std::accumulate(values_box.begin(), values_box.end(), 0.0) / values_box.size();
                vp.resize(s_box, titanlib::MV);
                for(int l = 0; l < s_box; l++) {
                    vp[l] = meanT;
                }
            }

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

            boost::numeric::ublas::vector<float> d(s_box);
            for(int i=0; i < s_box; i++) {
                d(i) = values_box[i] - vp[i]; // difference between actual temp and temperature from vertical profile
            }

            /* ---------------------------------------------------
               Beginning of real SCT
               ------------------------------------------------------*/
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
            for(int i=0; i<s_box; i++) {
                double acc = 0;
                for(int j=0; j<s_box; j++) {
                    acc += Sinv(i,j)*d(j);
                }
                Sinv_d(i) = acc;
            }
            for(int i=0; i<s_box; i++) {
                double acc = 0;
                for(int j=0; j<s_box; j++) {
                    acc += S(i,j)*Sinv_d(j);
                }
                ares_temp(i) = acc;
            }
            for(int i=0; i<s_box; i++) {
                Zinv(i) = (1/Sinv(i,i)); //Zinv<-1/diag(SRinv)
                ares(i) = ares_temp(i)-d(i); // ares<-crossprod(S,SRinv.d)-d[sel]
            }

            boost::numeric::ublas::vector<float> cvres(s_box);
            for(int i=0; i<s_box; i++) {
                cvres(i) = -1*Zinv(i) * Sinv_d(i);
            }

            double sig2o = 0;
            boost::numeric::ublas::vector<float> sig2o_temp(s_box), negAres_temp(s_box);
            for(int i=0; i<s_box; i++) {
                negAres_temp(i)=-1*ares(i);
                sig2o_temp(i) = d(i)*negAres_temp(i);
                sig2o += sig2o_temp(i);
            }

            sig2o = sig2o/s_box;
            if(sig2o < 0.01) {
                sig2o = 0.01;
            }

            // boost::numeric::ublas::vector<float> pog(s_box);
            // for(int i=0; i<s_box; i++) {
            //     pog(i) = cvres(i)*ares(i) / sig2o;
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
                    float pog = cvres(i) * ares(i) / sig2o;
                    assert(titanlib::is_valid(pog));
                    prob_gross_error[index] = std::max(pog, prob_gross_error[index]);
                    if((cvres(i) < 0 && pog > pos[index]) || (cvres(i) >= 0 && pog > neg[index])) {
                        flags[index] = 1;
                        thrown_out++;
                    }
                    checked[index] = 1;
                    ccount++;
                }
            }
            count_oi++;
        }
        if(thrown_out == 0) {
            if(iteration + 1 < num_iterations) {
                // std::cout << "Stopping early after " << iteration + 1<< " iterations" << std::endl;
            }
            break;
        }
        // std::cout << "Removing " << thrown_out << " Number of OI " << count_oi << std::endl;
        double e_time0 = titanlib::clock();
        // std::cout << e_time0 - s_time0 << std::endl;
    }

    return flags;
}
// end SCT //

//----------------------------------------------------------------------------//
// HELPER FUNCTIONS

void remove_flagged(ivec& indices, vec& distances, const ivec& flags) {
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
