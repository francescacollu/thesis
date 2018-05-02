//Corner-space renormalization method

#include <iostream>
#include <armadillo>
#include <iomanip>
#include <complex>
#include "init_matrices.hpp"

using namespace std;
using namespace arma;

int main() {
    
    double gamma = 1.;
    double Delta = 0.;
    int Niter = 4;
    int M = 10;
    
    //identity matrix
    cx_mat BlockI, I;
    
    //definition of the imaginary unit
    cx_double ii;
    
    //annihilation operator b's definition (fermionic case)
    cx_mat b, C;
    
    //definition of the hamiltonian
    cx_mat Sz, Sx, Sy, BlockH;
    
    init_matrices(Sx, Sy, Sz, I, b, ii, BlockH);
    
    //cout << "***Fin qui tutto bene1" << endl;
    
    cx_mat BlockSz = Sz;
    cx_mat BlockSx = Sx;
    cx_mat BlockSy = Sy;
    BlockI = I;
    C = b;
    cx_mat H = BlockH;
    
    cx_mat Hlink = kron(BlockSx, BlockSx) + kron(BlockSy, BlockSy);
    
    cx_vec Spectrum; //vettore caratterizzato dagli autovalori di L
    
    for(int i=1; i<Niter+1; i++)
    {
        //operatore liouvilliano: definito con una vettorizzazione con ordine per colonne (utilizzata da Armadillo)
        cx_mat first, second, third;
        cx_mat L;
        
        first = ii * trans(H) - gamma*0.5 * trans(C) * conj(C);
        second = -ii * H - gamma*0.5 * trans(C) * C;
        third = gamma * kron(conj(C), C);
        
        L = kron(first, BlockI) + kron(BlockI, second) + third;
        cout << "L = " << L.n_rows << " x " << L.n_cols << endl;
        
        //cx_vec Spectrum; //vettore costituito dagli autovalori di L
        cx_mat Evect; //matrice costituita dagli autovettori di L

        eig_gen(Spectrum, Evect, L);
        
        cx_vec rho_ss;
        cx_mat dm;
        
        for(int n=0; n<Spectrum.n_rows; n++)
        {
            if(Spectrum[n]==0.)
            {
                rho_ss = Evect.col(n);
            }
        }
        
        dm = reshape(rho_ss, sqrt(L.n_rows), sqrt(L.n_rows));
        cout << "Tr(DM) = " << trace (dm) << endl;
        
        
        //merging
        
        cx_vec spectrumDM, spectrumDmAUB;
        cx_mat eigenvecDM, eigenvecDmAUB;
        
        eig_gen(spectrumDM, eigenvecDM, dm);
        
        spectrumDmAUB = kron(spectrumDM, spectrumDM);
        eigenvecDmAUB = kron(eigenvecDM, eigenvecDM);
        
        int eigenvecCol = eigenvecDmAUB.n_cols;
        
        int Nkeep = min(eigenvecCol, M);
        cx_mat Omatr;
        Omatr = eigenvecDmAUB.tail_cols(Nkeep);
        
        BlockH = H;
        
        H = kron(BlockH, BlockI) + kron(BlockI, BlockH) + kron(BlockSx, BlockSx) + kron(BlockSy, BlockSy);
        BlockSx = kron(BlockI, BlockSx);
        BlockSy = kron(BlockI, BlockSy);
        C = kron(C, BlockI) + kron(BlockI, C);
        BlockI = kron(BlockI, BlockI);
        
        H = trans(Omatr)*H*Omatr;
        C = trans(Omatr)*C*Omatr;
        BlockSx = trans(Omatr)*BlockSx*Omatr;
        BlockSy = trans(Omatr)*BlockSy*Omatr;
        BlockI = trans(Omatr)*BlockI*Omatr;
        
        //cout << "H = " << H.n_rows << " x " << H.n_cols << endl;
        //cout << "C = " << C.n_rows << " x " << C.n_cols << endl;
        //cout << "Sx = " << BlockSx.n_rows << " x " << BlockSx.n_cols << endl;
        //cout << "Sy = " << BlockSy.n_rows << " x " << BlockSy.n_cols << endl;
        
        //cout << "***Fin qui tutto bene2" << endl;
        
        
    }
    
    cout << Spectrum << endl;
    
    
    return 0;
}
