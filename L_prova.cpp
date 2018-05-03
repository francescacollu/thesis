//prova di diagonalizzazione dell'operatore di Lindblad

#include <iostream>
#include <armadillo>
#include <iomanip>
#include <complex>

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
    Sx << cx_double(0., 0.) << cx_double(0.5, 0.) << endr << cx_double(0.5, 0.) << cx_double(0., 0.) << endr;
    Sy << cx_double(0., 0.) << cx_double(0., -0.5) << endr << cx_double(0., 0.5) << cx_double(0., 0.) << endr;
    
    cx_mat H, C, L, I, first, second, third;
    
    H = kron(Sx, kron(Sx, kron(BlockI, BlockI))) + kron(BlockI, kron(Sx, kron(Sx, BlockI))) + kron(BlockI, kron(BlockI, kron(Sx, Sx))) + kron(Sy, kron(Sy, kron(BlockI, BlockI))) + kron(BlockI, kron(Sy, kron(Sy, BlockI))) + kron(BlockI, kron(BlockI, kron(Sy, Sy)));
    //H = 0. * H;
    
    C = kron(b, kron(BlockI, kron(BlockI, BlockI))) + kron(BlockI, kron(b, kron(BlockI, BlockI))) + kron(BlockI, kron(BlockI, kron(b, BlockI))) + kron(BlockI, kron(BlockI, kron(BlockI, b)));
    
    I = kron(BlockI, kron(BlockI, kron(BlockI, BlockI)));
    
    first = ii * trans(H) - 0.5 * gamma * C.st() * C;
    second = - ii * H - 0.5 * gamma * trans(C) * conj(C);
    third = kron(C, conj(C));
    
    L = kron(first, I) + kron(I, second) + gamma * third;
    
    cx_vec Spectrum; //vettore costituito dagli autovalori di L
    //vec ReSpectrum, SortedReSpectrum;

    eig_gen(Spectrum, L);
    
    //ReSpectrum = real(Spectrum);
    //SortedReSpectrum = sort(ReSpectrum);
    
    cout << Spectrum << endl;
    
    double SpectrumRows = Spectrum.n_rows;
    
    cout << SpectrumRows << endl;
    
    
    
    
    
    
    /*cx_mat L, Hduo, Hquartet, Cduo, Cquartet;
    
    Hduo = kron(Sx, Sx) + kron(Sy, Sy); //4x4
    Hquartet = kron(kron(BlockI, BlockI), Hduo) + kron(Hduo, kron(BlockI, BlockI)) + kron(kron(BlockI, Sx), kron(Sx, BlockI)) + kron(kron(BlockI, Sy), kron(Sy, BlockI)); //16x16
    
    Cduo = gamma * (kron(C, BlockI) + kron(BlockI, C)); //4x4
    Cquartet = gamma * (kron(kron(BlockI, BlockI), Cduo) + kron(Cduo, kron(BlockI, BlockI))); //16x16
    
    cx_mat first, second, third;
    first = ii * trans(Hquartet) - gamma*0.5 * trans(Cquartet) * conj(Cquartet); //16x16
    first = kron(first, kron(kron(BlockI,BlockI), kron(BlockI, BlockI))); //16x16
    second = -ii * Hquartet - gamma*0.5 * trans(Cquartet) * Cquartet; //16x16
    second = kron((kron(kron(BlockI, BlockI), kron(BlockI, BlockI))), second); //16x16
    third = gamma * kron(conj(Cquartet), Cquartet); //256x256
    
    L = first + second + third; //256x256
    
    cx_vec Spectrum; //vettore costituito dagli autovalori di L
    
    eig_gen(Spectrum, L);
    
    cout << Spectrum << endl;*/
    
    return 0;
}
