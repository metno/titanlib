%module titanlib
%include "typemaps.i"
%include "std_vector.i"
%include "std_string.i"
namespace std {
  %template(IntVector) vector<int>;
  %template(DoubleVector) vector<double>;
  %template(FloatVector) vector<float>;
  %template(FloatVector2) vector<vector<float> >;
}
%apply std::vector<int>& OUTPUT { std::vector<int>& flags };
%apply std::vector<float>& OUTPUT { std::vector<float>& sct };
%apply std::vector<float>& OUTPUT { std::vector<float>& x_coords };
%apply std::vector<float>& OUTPUT { std::vector<float>& y_coords };
%apply std::vector<float>& OUTPUT { std::vector<float>& z_coords };
%apply std::vector<float>& OUTPUT { std::vector<float>& distances };
%{
/*  Put header files here or function declarations like below */
#include "titanlib.h"
%}
%include "titanlib.h"
