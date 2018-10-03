//prova di diagonalizzazione dell'operatore di Lindblad - 4 siti

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
    b << cx_double(0., 0.) << cx_double(0., 0.) << endr << cx_double(1., 0.) << cx_double(0., 0.) << endr;
    
    //Pauli matrices
    cx_mat Sx, Sy, Sz;
    Sx << cx_double(0., 0.) << cx_double(1., 0.) << endr << cx_double(1., 0.) << cx_double(0., 0.) << endr;
    Sy << cx_double(0., 0.) << cx_double(0., -1.) << endr << cx_double(0., 1.) << cx_double(0., 0.) << endr;
    Sz << cx_double(1., 0.) << cx_double(0., 0.) << endr << cx_double(0, 0.) << cx_double(-1., 0.) << endr;
    
    cx_mat H, C1, C2, C3, C4, L, I, first, second, third;
    
    H = - kron(Sx, kron(Sx, kron(BlockI, BlockI))) - kron(BlockI, kron(Sx, kron(Sx, BlockI))) - kron(BlockI, kron(BlockI, kron(Sx, Sx))) - 0.5 * (kron(Sy, kron(Sy, kron(BlockI, BlockI))) +  kron(BlockI, kron(Sy, kron(Sy, BlockI))) + kron(BlockI, kron(BlockI, kron(Sy, Sy))));
    
    C1 = kron(BlockI, kron(BlockI, kron(BlockI, b)));
    C2 = kron(BlockI, kron(BlockI, kron(BlockI, b)));
    C3 = kron(BlockI, kron(BlockI, kron(BlockI, b)));
    C4 = kron(BlockI, kron(BlockI, kron(BlockI, b)));
    
    I = kron(BlockI, kron(BlockI, kron(BlockI, BlockI)));
    
    first = ii * H.st() - 0.5 * gamma * C1.st() * conj(C1)- 0.5 * gamma * C2.st() * conj(C2)- 0.5 * gamma * C3.st() * conj(C3)- 0.5 * gamma * C4.st() * conj(C4);
    second = - ii * H - 0.5 * gamma * trans(C1) * C1 - 0.5 * gamma * trans(C2) * C2 - 0.5 * gamma * trans(C3) * C3- 0.5 * gamma * trans(C4) * C4;
    third = kron(conj(C1), C1) + kron(conj(C2), C2) + kron(conj(C3), C3) + kron(conj(C4), C4);
    
    L = kron(first, I) + kron(I, second) + gamma*third;
    
    cx_vec Spectrum; //vettore costituito dagli autovalori di L
    cx_mat Evect;
    
    eig_gen(Spectrum, Evect, L);
    
    cx_vec SortedSpectrum = sort(Spectrum);
    uvec indicesSpectrum = sort_index(Spectrum);
    
    cout << fixed << setprecision(10) << SortedSpectrum << endl;
    
    cx_mat EvectSorted;
    EvectSorted =  Evect.cols(indicesSpectrum);
    
    cx_vec rho_ss;
    cx_mat dm;
    
    int Lrows = L.n_rows;
    
    rho_ss = EvectSorted.col(0);
    dm = reshape(rho_ss, sqrt(Lrows), sqrt(Lrows));
    
    cout << "La traccia di dm è " << trace(dm) << endl;
    
    dm = dm/trace(dm);
    cout << "La traccia di dm (giusta) è " << trace(dm) << endl;
    
    cx_mat SigmaZ;
    SigmaZ = kron(kron(Sz, BlockI), kron(BlockI, BlockI));
    cout << "Tr(p*SigmaZ) = " << trace(dm*SigmaZ) << endl;
    return 0;
}
