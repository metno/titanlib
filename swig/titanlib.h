#ifndef TITANLIB_H
#define TITANLIB_H
#include <iostream>
#include <vector>
typedef std::vector<float> vec;

float calc_gamma(float shape, float scale);
bool sct(const vec lats,
        const vec lons,
        const vec elevs,
        int nmin,
        int nmax,
        int nminprof,
        float dzmin,
        float dhmin,
        float dz,
        const vec t2pos,
        const vec t2neg,
        const vec eps2,
        vec& flags);

bool plausibility(const vec lats,
        const vec lons,
        const vec elevs,
        int unixtime,
        float min,
        float max,
        vec& flags);

#endif
