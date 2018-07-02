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
    int Niter = 3;
    int M = 17;
    
    //identity matrix
    cx_mat BlockI, I;
    
    //definition of the imaginary unit
    cx_double ii;
    
    //annihilation operator b's definition (fermionic case)
    cx_mat b, C, BlockC;
    
    //definition of the hamiltonian
    cx_mat Sz, Sx, Sy, BlockH;
    
    init_matrices(Sx, Sy, Sz, I, b, ii, BlockH);
    
    cx_mat BlockSz = Sz;
    cx_mat BlockSxsx = Sx;
    cx_mat BlockSxdx = Sx;
    cx_mat BlockSysx = Sy;
    cx_mat BlockSydx = Sy;
    BlockI = I;
    C = b;
    cx_mat H = BlockH;
    
    cx_vec Spectrum, SpectrumSorted; //vettore degli autovalori di L
    cx_mat dm;
    
        
    for(int i=1; i<=Niter; i++)
    {
        //operatore liouvilliano: definito con una vettorizzazione con ordine per colonne (utilizzata da Armadillo)
        cx_mat first, second, third;
        cx_mat L;
        
        /*first = -ii * H - gamma*0.5 * trans(C) * C;
        second = ii * trans(H) - gamma*0.5 * trans(C) * conj(C);
        third = gamma * kron(C, conj(C));*/
        //L'espressione del liouvilliano scritta sopra è quella ottenuta effettuando una vettorizzazione per righe; poiché Armadillo fa una vettorizzazione per colonne, usiamo l'espressione scritta di seguito. N.B. Gli autovalori restano uguali, probabilmente per via della scelta di C => kron(conj(C), C) = kron(C, conj(C)).
        
        first = cx_double(0., 1.) * BlockH.st() - 0.5 * (C.st() * conj(C));
        first = kron(first, BlockI);
        second = cx_double(0., -1.) * BlockH - 0.5 * (trans(C) * C);
        second = kron(BlockI, second);
        third = kron(conj(C), C);
        
        L = first + second + third;
        cout << "L = " << L.n_rows << " x " << L.n_cols << endl;
        
        cx_mat Evect; //matrice costituita dagli autovettori di L
        
        eig_gen(Spectrum, Evect, L);
        
        SpectrumSorted = sort(Spectrum);
        uvec indicesSpectrum = sort_index(Spectrum);
        
        cx_mat EvectSorted;
        EvectSorted =  Evect.cols(indicesSpectrum);
        
        cx_vec rho_ss;
        
        //Testato che non ha rilevanza sul calcolo degli autovalori di L successivi.
        /*for(int n=0; n<Spectrum.n_rows; n++)
        {
            if(std::norm(Spectrum[n])<1.E-10)
            {
                Spectrum[n] = 0.;
            }
        }*/
        
        cout << "Spectrum di L: " << endl << fixed << setprecision(9) << SpectrumSorted << endl;
        int Lrows = L.n_rows;

        rho_ss = EvectSorted.col(0);
        dm = reshape(rho_ss, sqrt(Lrows), sqrt(Lrows));
        
        //merging
        cx_vec spectrumDM, spectrumDmAUB, phiphi, spectrumDmAUBsorted;
        cx_mat eigenvecDM, eigenvecDmAUB, eigenvecDmAUBsorted;
        
        eig_gen(spectrumDM, eigenvecDM, dm);
        
        spectrumDmAUB = kron(spectrumDM, spectrumDM);
        eigenvecDmAUB = kron(eigenvecDM, eigenvecDM);
        
        /*int k=0;
        for(int i=0; i<eigenvecDM.n_cols; i++)
        {
            for(int j=0; j<eigenvecDM.n_cols; j++)
            {
                cx_vec eveci = eigenvecDM.col(i);
                cx_vec evecj = eigenvecDM.col(j);
                
                phiphi = kron(eveci, evecj);
                eigenvecDmAUB.insert_cols(k, phiphi);
                k=k+1;
            }
        }*/
        
        spectrumDmAUBsorted = sort(spectrumDmAUB);
        uvec indices = sort_index(spectrumDmAUB);
        
        eigenvecDmAUBsorted =  eigenvecDmAUB.cols(indices);
        
        int eigenvecCol = eigenvecDmAUBsorted.n_cols;
        vec choice;
        choice << eigenvecCol << M << endr;
        
        int Nkeep = min(choice); //valore minimo del vettore choice
        cx_mat Omatr;
        Omatr = eigenvecDmAUBsorted.tail_cols(Nkeep);
        
        cout << "Omatr dim: " << Omatr.n_rows << " x " << Omatr.n_cols << endl;
        
        // rinormalizzazione
        H = kron(BlockH, BlockI) + kron(BlockI, BlockH) - kron(BlockSxsx, BlockSxdx) - cx_double(0.5, 0.) * kron(BlockSysx, BlockSydx);
        
        BlockSxsx = kron(BlockI, BlockSxsx);
        BlockSxdx = kron(BlockSxdx, BlockI);
        BlockSysx = kron(BlockI, BlockSysx);
        BlockSydx = kron(BlockSydx, BlockI);
        C = kron(C, BlockI) + kron(BlockI, C);
        BlockI = kron(BlockI, BlockI);
        
        BlockH = trans(Omatr)*H*Omatr;
        C = trans(Omatr)*C*Omatr;
        BlockSxsx = trans(Omatr)*BlockSxsx*Omatr;
        BlockSxdx = trans(Omatr)*BlockSxdx*Omatr;
        BlockSysx = trans(Omatr)*BlockSysx*Omatr;
        BlockSydx = trans(Omatr)*BlockSydx*Omatr;
        BlockI = trans(Omatr)*BlockI*Omatr;
    }
    
    return 0;
}
