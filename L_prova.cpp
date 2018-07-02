//prova di diagonalizzazione dell'operatore di Lindblad

#include <iostream>
#include <armadillo>
#include <iomanip>
#include <complex>
#include <fstream>

using namespace std;
using namespace arma;

int main() {
    
    double gamma = 1.;
    
    //identity matrix
    cx_mat BlockI;
    BlockI << cx_double(1., 0.) << cx_double(0., 0.) << endr << cx_double(0., 0.) << cx_double(1., 0.) << endr;
    
    //definition of the imaginary unit
    cx_double ii = cx_double(0., 1.);
    
    //annihilation operator
    cx_mat b;
    b << cx_double(0., 0.) << cx_double(0., 0.) << endr << cx_double(2., 0.) << cx_double(0., 0.) << endr;
    
    //Pauli matrices
    cx_mat Sx, Sy;
    Sx << cx_double(0., 0.) << cx_double(1., 0.) << endr << cx_double(1., 0.) << cx_double(0., 0.) << endr;
    Sy << cx_double(0., 0.) << cx_double(0., -1.) << endr << cx_double(0., 1.) << cx_double(0., 0.) << endr;
    
    cx_mat H, C, L, I, first, second, third;
    
    H = - kron(Sx, kron(Sx, kron(BlockI, BlockI))) - kron(BlockI, kron(Sx, kron(Sx, BlockI))) - kron(BlockI, kron(BlockI, kron(Sx, Sx))) -  0.5*(kron(Sy, kron(Sy, kron(BlockI, BlockI))) +  kron(BlockI, kron(Sy, kron(Sy, BlockI))) + kron(BlockI, kron(BlockI, kron(Sy, Sy))));
    
    C = kron(b, kron(BlockI, kron(BlockI, BlockI))) + kron(BlockI, kron(b, kron(BlockI, BlockI))) + kron(BlockI, kron(BlockI, kron(b, BlockI))) + kron(BlockI, kron(BlockI, kron(BlockI, b)));
    
    I = kron(BlockI, kron(BlockI, kron(BlockI, BlockI)));
    
    first = ii * trans(H) - 0.5 * gamma * C.st() * conj(C);
    second = - ii * H - 0.5 * gamma * trans(C) * C;
    third = kron(conj(C), C);
    
    L = kron(first, I) + kron(I, second) + gamma * third;
    
    cx_vec Spectrum; //vettore costituito dagli autovalori di L
    //vec ReSpectrum, SortedReSpectrum;

    eig_gen(Spectrum, L);
    
    //ReSpectrum = real(Spectrum);
    //SortedReSpectrum = sort(ReSpectrum);
    
    Spectrum = sort(Spectrum);
    
    
    double SpectrumRows = Spectrum.n_rows;
    ofstream myfile ("L_prova.dat");
    
    for(int j=0; j<SpectrumRows; j++)
    {
        myfile << fixed << std::scientific << real(Spectrum[j]) << "\t" << imag(Spectrum[j]) << endl;
    }
    
    return 0;
}
