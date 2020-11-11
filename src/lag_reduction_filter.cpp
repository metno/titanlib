#include "titanlib.h"

using namespace titanlib;

vec titanlib::lag_reduction_filter(const vec& times, const vec& values, float a, float b, float k1, float k2, int n) {
    vec result(times.size());
    result[0] = values[0];
    for(int i = 1; i < times.size(); i++) {
        float curr = result[i - 1];
        float f0 = values[i - 1];
        float f1 = values[i];
        float prev = curr;
        float time_diff = times[i] - times[i - 1];
        float deriv1 = (f1 - f0) / time_diff;
        float deriv2 = 0; // # (f1 - f0) - (f0 - fm1);
        float dt = 1.0 / n;
        for(int ti = 0; ti < n; ti++) {
            float t = (ti + 1) * dt;
            float fi = f0 + t * (f1 - f0);
            float F0 = 1.0 / (a*k1 + b*k2) * (deriv2 + (k1 + k2) * deriv1 + k1 * k2 * (fi - curr));
            fi = f0 + (t+dt) * (f1 - f0);
            float F1 = 1.0 / (a*k1 + b*k2) * (deriv2 + (k1 + k2) * deriv1 + k1 * k2 * (fi - curr));
            curr = prev + dt / 2 * (F0 + F1);
            prev = curr;
        }
        result[i] = curr;
    }
    return result;
}

