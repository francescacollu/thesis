#pragma once
#include <armadillo>

namespace corner {
    void check(const bool& condition, const char* class_name, const char* msg);
    void info(const bool& condition, const char* class_name, const char* msg);
    void warning(const bool& condition, const char* class_name, const char* msg);
    bool approx_equal(arma::cx_double, arma::cx_double, double tol=1E-10);
}
