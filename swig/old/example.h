#ifndef EXAMPLE_H
#define EXAMPLE_H
#include <vector>
      enum variable { temperature, precipitation};
      int fact(int n);
      int my_mod(int x, int y);
      char *get_time();
      float get_gamma();
      float get_gamma(float shape);
      float test(int a, int b, int c=3);
      std::vector<float> get_mean(const std::vector<float> val, std::vector<float>& output);
#endif
