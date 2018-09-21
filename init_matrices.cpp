#include "init_matrices.hpp"

using namespace arma;

void init_matrices (arma::cx_mat& Sx, arma::cx_mat& Sy, arma::cx_mat& Sz, arma::cx_mat& I, arma::cx_mat& b, arma::cx_double& ii, arma::cx_mat& BlockH)
{
    I << cx_double(1., 0.) << cx_double(0., 0.) << endr << cx_double(0., 0.) << cx_double(1., 0.) << endr;
    
    ii = cx_double(0., 1.);
    
    b << cx_double(0., 0.) << cx_double(0., 0.) << endr << cx_double(1., 0.) << cx_double(0., 0.) << endr;
    
    Sz << cx_double(1., 0.) << cx_double(0., 0.) << endr << cx_double(0., 0.) << cx_double(-1., 0.) << endr;
    Sx << cx_double(0., 0.) << cx_double(1., 0.) << endr << cx_double(1., 0.) << cx_double(0., 0.) << endr;
    Sy << cx_double(0., 0.) << cx_double(0., -1.) << endr << cx_double(0., 1.) << cx_double(0., 0.) << endr;
    
    BlockH << cx_double(0., 0.) << cx_double(0., 0.) << endr << cx_double(0., 0.) << cx_double(0., 0.) << endr;
}
