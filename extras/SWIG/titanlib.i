%module titanlib
%init %{
#if defined(SWIGPYTHON)
    import_array();
#endif
    titanlib::initialize_omp();
%}
%include exception.i
/*
      SWIG_MemoryError
      SWIG_IOError
      SWIG_RuntimeError
      SWIG_IndexError
      SWIG_TypeError
      SWIG_DivisionByZero
      SWIG_OverflowError
      SWIG_SyntaxError
      SWIG_ValueError
      SWIG_SystemError
*/
%exception {
    try {
        $action
    }
    catch (std::invalid_argument &e) {
        std::string s(e.what());
        SWIG_exception(SWIG_ValueError, s.c_str());
    }
    catch (titanlib::not_implemented_exception &e) {
        std::string s(e.what());
        SWIG_exception(SWIG_RuntimeError, s.c_str());
    }
    catch (std::exception &e) {
        std::string s(e.what());
        SWIG_exception(SWIG_RuntimeError, s.c_str());
    }
    catch (...) {
        SWIG_exception(SWIG_RuntimeError, "Unknown exception");
    }
}
%include "vector.i"
%apply std::vector<int>& OUTPUT { std::vector<int>& flags };
%apply std::vector<float>& OUTPUT { std::vector<float>& prob_gross_error };
%apply std::vector<float>& OUTPUT { std::vector<float>& rep };
%apply std::vector<int>& OUTPUT { std::vector<int>& boxids };
%apply std::vector<float>& OUTPUT { std::vector<float>& x_coords };
%apply std::vector<float>& OUTPUT { std::vector<float>& y_coords };
%apply std::vector<float>& OUTPUT { std::vector<float>& z_coords };
%apply std::vector<float>& OUTPUT { std::vector<float>& distances };
%apply std::vector<std::vector<float> >& OUTPUT { std::vector<std::vector<float> >& distances };
%apply std::vector<float>& OUTPUT { std::vector<float>& scores };
%apply std::vector<float>& OUTPUT { std::vector<float>& num_inner };
%apply std::vector<float>& OUTPUT { std::vector<float>& horizontal_scale };
%apply int& OUTPUT { int& X1_out };
%apply int& OUTPUT { int& Y1_out };
%apply int& OUTPUT { int& X2_out };
%apply int& OUTPUT { int& Y2_out };

%{
#include "titanlib.h"
%}

%include "titanlib.h"
