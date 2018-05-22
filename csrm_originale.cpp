//Corner-space renormalization method

#include <iostream>
#include <armadillo>
#include <iomanip>
#include <complex>
#include "init_matrices.hpp"

using namespace std;
using namespace arma;

//definisco la classe eigen; ogni eigen è definito da un autovalore e un autovettore.
/*class eigen{
public:
    cx_double value;
    cx_vec vector;
    
    eigen() {value = 0;};
    eigen (cx_double, cx_vec);
    eigen operator* (const eigen&);
};

eigen::eigen (cx_double a, cx_vec b){
    value = a;
    vector = b;
}

eigen eigen::operator* (const eigen& param)
{
    eigen temp;
    temp.value = value*param.value;
    temp.vector = kron(vector, param.vector);
    return(temp);
}*/

int main() {
    
    double gamma = 1.;
    double Delta = 0.;
    int Niter = 3;
    //int M = 4;
    
    //identity matrix
    cx_mat BlockI, I;
    
    //definition of the imaginary unit
    cx_double ii;
    
    //annihilation operator b's definition (fermionic case)
    cx_mat b, C;
    
    //definition of the hamiltonian
    cx_mat Sz, Sx, Sy, BlockH, sigmaZ;
    
    init_matrices(Sx, Sy, Sz, I, b, ii, BlockH, sigmaZ);
    
    //cout << "***Fin qui tutto bene1" << endl;
    
    cx_mat BlockSz = Sz;
    cx_mat BlockSxsx = Sx;
    cx_mat BlockSxdx = Sx;
    cx_mat BlockSysx = Sy;
    cx_mat BlockSydx = Sy;
    BlockI = I;
    C = b;
    cx_mat H = BlockH;
    
    cx_vec Spectrum; //vettore caratterizzato dagli autovalori di L
    cx_mat dm, dmY;
    
    ofstream myfile("sigmaZ.dat");
    
    for(int M=4; M<25; M++)
    {
        
    for(int i=1; i<Niter+1; i++)
    {
        //operatore liouvilliano: definito con una vettorizzazione con ordine per colonne (utilizzata da Armadillo)
        cx_mat first, second, third;
        cx_mat L;
        
        /*first = -ii * H - gamma*0.5 * trans(C) * C;
        second = ii * trans(H) - gamma*0.5 * trans(C) * conj(C);
        third = gamma * kron(C, conj(C));*/
        //L'espressione del liouvilliano scritta sopra è quella ottenuta effettuando una vettorizzazione per righe; poiché Armadillo fa una vettorizzazione per colonne, usiamo l'espressione scritta di seguito. N.B. Gli autovalori restano uguali, probabilmente per via della scelta di C => kron(conj(C), C) = kron(C, conj(C)).
        
        first = ii * trans(H) - gamma*0.5 * C.st() * conj(C);
        second = -ii * H - gamma*0.5 * trans(C) * C;
        third = gamma * kron(conj(C), C);
        
        L = kron(first, BlockI) + kron(BlockI, second) + third;
        cout << "L = " << L.n_rows << " x " << L.n_cols << endl;
        
        //cout << "L: " << endl << L << endl;
        
        cx_mat Evect; //matrice costituita dagli autovettori di L
        
        eig_gen(Spectrum, Evect, L);
        
        //cout << "Spectrum di L: " << endl << fixed << setprecision(9) << Spectrum << endl;
        //cout << Evect << endl;
        
        cx_vec rho_ss;
        //cx_mat dm;
        
        /*for(int n=0; n<Spectrum.n_rows; n++)
        {
            if(std::norm(Spectrum[n])==0. || std::norm(Spectrum[n])<1.E-10)
            {
                //cout << "n = " << n << endl;
                //cout << Spectrum[n] << endl;
                rho_ss = Evect.col(n);
                dm = reshape(rho_ss, sqrt(L.n_rows), sqrt(L.n_rows));
                if(norm(trace(dm))==1. || (norm(trace(dm))<1.1 && norm(trace(dm))>0.9))
                {
                    dmY = dm;
                    cout << "Trace of the density matrix = " << trace(dmY) << endl;
                }
            }
        }*/
        
        int Lrows = L.n_rows;

        rho_ss = Evect.col(0);
        dm = resize(rho_ss, sqrt(Lrows), sqrt(Lrows)); //ho cambiato il reshape in resize!!!
        
        //cout << dm << endl;
        
        //cout << rho_ss << endl << endl;
        
        //cout << "<SigmaZ> = " << trace(dm*sigmaZ) << endl;
        
        //merging
        cx_vec spectrumDM, spectrumDmAUB, phiphi, spectrumDmAUBsorted;
        cx_mat eigenvecDM, eigenvecDmAUB, eigenvecDmAUBsorted;
        
        eig_gen(spectrumDM, eigenvecDM, dm);
        
        spectrumDmAUB = kron(spectrumDM, spectrumDM);
        
        int k=0;
        
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
        }
        
        //cout << "matrice del kron prod dei phi:" << endl << eigenvecDmAUB << endl;
        
        //cout << "KronProd degli eigenvecDM:" << endl << eigenvecDmAUB << endl;
        
        //cout << eigenvecDmAUB << endl;
        
        spectrumDmAUBsorted = sort(spectrumDmAUB);
        uvec indices = sort_index(spectrumDmAUB);
        
        eigenvecDmAUBsorted =  eigenvecDmAUB.cols(indices);
        
        //cout << spectrumDmAUBsorted << endl;
        //cout << eigenvecDmAUB << endl;
        
        int eigenvecCol = eigenvecDmAUBsorted.n_cols;
        
        vec choice;
        choice << eigenvecCol << M << endr;
        
        int Nkeep = min(choice);
        
        //cout << "Nkeep = " << Nkeep << endl;
        cx_mat Omatr;
        Omatr = eigenvecDmAUBsorted.tail_cols(Nkeep);
        
        cout << "Omatr dim: " << Omatr.n_rows << " x " << Omatr.n_cols << endl;
        
        BlockH = H;
        
        H = kron(BlockH, BlockI) + kron(BlockI, BlockH) - kron(BlockSxsx, BlockSxdx) - 0.5 * kron(BlockSysx, BlockSydx);
        
        BlockSxsx = kron(BlockI, BlockSxsx);
        BlockSxdx = kron(BlockSxdx, BlockI);
        BlockSysx = kron(BlockI, BlockSysx);
        BlockSydx = kron(BlockSydx, BlockI);
        sigmaZ = kron(BlockI, sigmaZ);
        C = kron(C, BlockI) + kron(BlockI, C);
        BlockI = kron(BlockI, BlockI);
        
        H = trans(Omatr)*H*Omatr;
        C = trans(Omatr)*C*Omatr;
        BlockSxsx = trans(Omatr)*BlockSxsx*Omatr;
        BlockSxdx = trans(Omatr)*BlockSxdx*Omatr;
        BlockSysx = trans(Omatr)*BlockSysx*Omatr;
        BlockSydx = trans(Omatr)*BlockSydx*Omatr;
        BlockI = trans(Omatr)*BlockI*Omatr;
        sigmaZ = trans(Omatr)*sigmaZ*Omatr;
        
        //cout << "***Fin qui tutto bene2" << endl;
        
    }
    
        myfile << trace(dm*sigmaZ) << endl;
    
    //cout << "sigmaZ dimensions: " << sigmaZ.n_rows << " x " << sigmaZ.n_cols << endl;
    //cout << "dm dimensions: " << dm.n_rows << " x " << dm.n_cols << endl;
    //cout << "<SigmaZ> = " << fixed << setprecision(9) << trace(dm*sigmaZ) << endl;
    //cout << fixed << setprecision(30) << Spectrum << endl;
    }
    
    
    /*int SpectrumRows = Spectrum.n_rows;
    ofstream myfile("CSRM.dat");
    
    for(int j=0; j<SpectrumRows; j++)
    {
        myfile << fixed << setprecision(20) << real(Spectrum[j]) << "\t" << imag(Spectrum[j]) << endl;
    }*/
    
    return 0;
    
}
