/*  File : example.c */

#include <time.h>
#include <boost/math/distributions/gamma.hpp>
#include "example.h"
#include <iostream>
int fact(int n) {
    if (n <= 1) return 1;
    else return n*fact(n-1);
}

int my_mod(int x, int y) {
    return (x%y);
}

std::vector<float> get_mean(const std::vector<float> val, std::vector<float>& output) {
    float count = 0;
    std::cout << "Number of values: " << val.size() << std::endl;
    output.resize(val.size());
    for(int i = 0; i < val.size(); i++) {
        count += val[i];
        output[i] = val[i] * 2;
    }
    return val;
}

float get_gamma() {
    float shape = 1;
    float scale = 2;
    boost::math::gamma_distribution<> dist(shape, scale);
    float value = boost::math::quantile(dist, 0.5);
    return value;
}

float get_gamma(float shape) {
    float scale = 2;
    boost::math::gamma_distribution<> dist(shape, scale);
    float value = boost::math::quantile(dist, 0.5);
    return value;
}
float test(int a, int b, int c) {
    std::cout << a << " " << b << " " << c << std::endl;
    return 1;
}
char *get_time()
{
    time_t ltime;
    time(&ltime);
    return ctime(&ltime);
}
