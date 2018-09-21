/*
Class describing an hermitian matrix
*/

#pragma once

#include<armadillo>

class HMatrix
{

	arma::cx_mat m;
	arma::cx_vec eigval;
	arma::cx_mat eigvec;

    arma::cx_vec Projection(arma::cx_vec oldvec, arma::cx_vec newvec);
public:
	HMatrix(const arma::cx_mat&);
	int size() const;
	// Build an orthonormal basis from the hermitian matrix
	arma::cx_mat GetONBasis();

	// Get order of degeneration of i-th eigenvalue
	int GetDegeneration(const int& i);

	// Apply Gram-Schmidt algorithm
    arma::cx_mat GramSchmidt(const arma::cx_mat&);

    //Get vector of sorted eigenvalues
    const arma::cx_vec& GetEigenvalues() const;
};

