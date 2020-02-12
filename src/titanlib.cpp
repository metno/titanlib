#include <iostream>
#include <vector>
#include "titanlib.h"
#include <boost/math/distributions/gamma.hpp>

float titanlib::calc_gamma(float shape, float scale) {
    boost::math::gamma_distribution<> dist(shape, scale);
    float value = boost::math::quantile(dist, 0.5);
    return value;
}
