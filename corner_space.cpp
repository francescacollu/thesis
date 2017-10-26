//Corner Space Renormalization Method (CSRM)
//L = (-i h - gamma * C+ * C) kron 1 + 1 kron (i h+ - gamma * C+ * C*) + gamma (C kron C*)

#include <iostream>
#include <armadillo>
#include <iomanip>
#include <complex>

using namespace std;
using namespace arma;

int main() {
    
    int m = 10;
    int Niter = 10;
    double gamma = 0.5;
    double Delta = 0.;
    int n = 4; //number of sites
    
    mat I = eye<mat>(2,2);
    mat Ifourfour = eye<mat>(4,4);
    
    //definition of the imaginary unit
    complex <double> ii (0, 1.);
    
    //annihilation operator's b definition
    mat b = zeros<mat>(n, n);
    for (int i=1; i < n; i++)
    {
        b(i-1, i) = sqrt(i);
    }
    
    //jumps operator's C definition
    mat C = gamma * b;
    
    
    //definition of the hamiltonian
    mat Sz, Sp, Sm;
    Sz << 0.5 << 0 << endr << 0 << -0.5 << endr;
    Sp << 0 << 0 << endr << 1 << 0 << endr;
    Sm << 0 << 1 << endr << 0 << 0 << endr;

    mat BlockSz = Sz;
    mat BlockSp = Sp;
    mat BlockSm = Sm;
    mat BlockI = I;
    mat BlockH = zeros<mat>(2,2);
    mat H_super;
    double Energy = 0.;
    double LastEnergy;
    
    //cout << "Sistem Size" << "\t" << "Energy" << "\t" << "Energy per bond" << "\t" << "Ener2" << endl;
    
    
    for (int i=1; i<Niter+1; i++)
    {
        
        double SystSize, Ener2, EnergyPerBond;
        SystSize = pow(2, 2*i);
        
        H_super = kron(BlockH, BlockI) + kron(BlockI, BlockH) - Delta * kron(BlockSz, BlockSz) + 0.5 * (kron(BlockSp, BlockSm) + kron(BlockSm, BlockSp));
        H_super = 0.5 * (H_super + trans(H_super));
        
        //definition of the Lindblad operator
        cx_mat first = -ii * H_super - gamma * trans(C) * C;
        //cout << kron(first, Ifourfour) << endl;
        
        cx_mat second = ii * H_super - gamma * trans(C) * C;
        //cout << kron(Ifourfour, second) << endl;
        
        mat third = gamma * kron(C, conj(C));
        //cout << third << endl;
        
        cx_mat L = kron(first, Ifourfour) + kron(Ifourfour, second) + third;
        
        L = 0.5 * (L + trans(L));
        
        LastEnergy = Energy;
        
        //diagonalization
        vec Spectrum;
        cx_mat Evect;
        eig_sym(Spectrum, Evect, L);
        
        Energy = Spectrum[0];
        Ener2 = Energy/SystSize;
        EnergyPerBond = (Energy-LastEnergy)/(SystSize/2);
        
        
        //truncation
        int Nkeep, SpectrumRows;
        cx_mat Omatr;
        
        SpectrumRows = Spectrum.n_rows;
        
        Nkeep = min(SpectrumRows, m);
        
        Omatr = Evect.cols(0, Nkeep-1);
        
        cout << fixed << setprecision(6) << SystSize << "\t" << Energy << "\t" << EnergyPerBond << "\t" << Ener2 << endl;
        
        /*BlockSz = kron(BlockI, BlockSz);
        BlockSp = kron(BlockI, BlockSp);
        BlockSm = kron(BlockI, BlockSm);
        BlockI = kron(BlockI, BlockI);
        
        BlockH = trans(Omatr)*H_super*Omatr;
        BlockSz = trans(Omatr)*BlockSz*Omatr;
        BlockSp = trans(Omatr)*BlockSp*Omatr;
        BlockSm = trans(Omatr)*BlockSm*Omatr;
        BlockI = trans(Omatr)*BlockI*Omatr;
        L = trans(Omatr)*L*Omatr;
         */
    }
    return 0;
}
