// corner-space-renormalization-method
// codice strutturato in classi
// 2 Lindblad reservoirs in the 1st and in the 8th sites.

#include <iostream>
#include <armadillo>
#include <iomanip>
#include <complex>
#include <vector>
#include <fstream>
#include "init_matrices_prova.hpp"
#include "HMatrix_prova.hpp"
#include "utility_prova.hpp"
#include "LiouvSO_prova.hpp"

using namespace std;
using namespace arma;
using namespace corner;

int main()
{
    double gamma = 1.;
    double Jx = 1.;
    double Jy = 0.5;
    double Jz = 1.;
    int Niter = 3;
    int M = 40;
    ofstream myfile("expValueSz_2reservoir.txt");
    
    cx_mat BlockI, I; //identity matrix
    cx_double ii; //imaginary unit
    cx_mat b; //jump operator
    cx_mat Sx, Sy, Sz, BlockH; //H's matrices
    init_matrices(Sx, Sy, Sz, I, b, ii, BlockH);
    cx_mat BlockSxsx = Sx;
    cx_mat BlockSxdx = Sx;
    cx_mat BlockSysx = Sy;
    cx_mat BlockSydx = Sy;
    cx_mat BlockSzsx = Sz;
    cx_mat BlockSzdx = Sz;
    cx_mat H = BlockH;
    BlockI = I;
    
    vector<cx_mat> sigmaZ(pow(2, Niter));
    vector<cx_mat> sigmaZnew(pow(2, Niter));
    vector<cx_mat> sigmaZnewnew(pow(2, Niter));
    
    sigmaZ[0] = Sz;
    
    cx_mat Da, Db, D;
    vector<cx_mat> Ca(pow(2,Niter)/2), Cb(pow(2,Niter)/2), C(pow(2,Niter));
    vector<cx_mat> CnewA(pow(2,Niter)/2), CnewB(pow(2,Niter)/2), Cnew(pow(2,Niter));
    vector<cx_mat> CnewnewA(pow(2,Niter)/2), CnewnewB(pow(2,Niter)/2), Cnewnew(pow(2,Niter));
    vector<cx_mat> DsliceA(pow(2,Niter)/2), DsliceB(pow(2,Niter)/2), Dslice(pow(2,Niter));
    
    Ca[0] = b;
    Cb[0] = trans(b);
    
    Da = -0.5*gamma*kron(Ca[0].st()*conj(Ca[0]),BlockI) - 0.5*gamma*kron(BlockI,trans(Ca[0])*Ca[0]) + kron(conj(Ca[0]),Ca[0]);
    Db = -0.5*gamma*kron(Cb[0].st()*conj(Cb[0]),BlockI) - 0.5*gamma*kron(BlockI,trans(Cb[0])*Cb[0]) + kron(conj(Cb[0]),Cb[0]);
    
    cx_mat L_A, L_B, L_AUB; //liouvillian superoperator
    
    for(int i=0; i<Niter; i++)
    {
        ////////// First block //////////
        L_A = ii*kron(BlockH.st(), BlockI) - ii*kron(BlockI, BlockH) + Da;
        //cout << "L_Adim: " << L_A.n_rows << "x" << L_A.n_cols << endl;
        cout << "***DM A***" << endl;
        LiouvSO LA(L_A);
        cx_mat dmA = LA.GetDM();
        HMatrix DM_A(dmA);
        cx_mat EvectDmA = DM_A.GetONBasis();
        cx_vec SpectrumDmA = DM_A.GetEigenvalues();
        
        
        ////////// Second block //////////
        L_B = ii*kron(BlockH.st(), BlockI) - ii*kron(BlockI, BlockH) + Db;
        //cout << "L_Bdim: " << L_B.n_rows << "x" << L_B.n_cols << endl;
        cout << "***DM B***" << endl;
        LiouvSO LB(L_B);
        cx_mat dmB = LB.GetDM();
        HMatrix DM_B(dmB);
        cx_mat EvectDmB = DM_B.GetONBasis();
        cx_vec SpectrumDmB = DM_B.GetEigenvalues();
        
        //*//*//*//*// A U B //*//*//*//*//
        ////////// Building the pseudo-orthogonal matrix/////////
        cx_mat dmAUB = kron(dmA, dmB);
        HMatrix DM_AUB(dmAUB);
        cx_mat EvectDmAUB = DM_AUB.GetONBasis();
        cx_vec SpectrumDmAUB = DM_AUB.GetEigenvalues();
        cout << "dmRows = " << dmAUB.n_rows << "x" << dmAUB.n_cols << endl;
        
        int EvectCols = EvectDmAUB.n_cols;
        vec choice;
        choice << EvectCols << M << endr;
        int Nkeep = min(choice); //valore minimo del vettore choice
        
        cx_mat O = EvectDmAUB.tail_cols(Nkeep); //matrice pseudo-ortogonale
        
        
        //*//*//*//*//* MERGE *//*//*//*//*//
        ////////// Merging the Hamiltonians //////////
        cout << "Merging the hamiltonians..." << endl;
        H = kron(BlockH, BlockI) + kron(BlockI, BlockH) - Jx*kron(BlockSxsx, BlockSxdx) - Jy*kron(BlockSysx, BlockSydx) - Jz*kron(BlockSzsx, BlockSzdx);
        
        cout << "Working on the dissipators..." << endl;
        // only the first block has (one) dissipator (in the 0 site)
        for(int k=0; k<pow(2, i+1)/2; k++) // distribution of the dissipator(s) over the first block-for the next loop
        {
            if(k==0)
            {
                CnewA[k]=kron(Ca[k], BlockI);
            }
            else
            {
                CnewA[k] = kron(BlockI, BlockI); //Ca[k-pow(2,i)+1]);
            }
        }
        
        cout << "Working on the dissipators..." << endl;
        for(int k=0; k<pow(2, i+1)/2; k++) // distribution of the dissipator(s) over the second block-for the next loop
        {
            if(k==(pow(2, i+1)/2)-1)
            {
                CnewB[k]=kron(BlockI, Cb[k-pow(2,i)+1]);
            }
            else
            {
                CnewB[k] = kron(BlockI, BlockI);//C[k-pow(2,i)]);
            }
        }
        
        cout << "Renormalization..." << endl;
        ////////// Renormalization //////////
        for(int k=0; k<pow(2, i+1); k++)
        {
            if(k<pow(2,i))
            {
                sigmaZnew[k]=kron(sigmaZ[k], BlockI);
            }
            else
            {
                sigmaZnew[k] = kron(BlockI, sigmaZ[k-pow(2,i)]);
            }
        }
        
        BlockSxsx = kron(BlockI, BlockSxsx);
        BlockSxdx = kron(BlockSxdx, BlockI);
        BlockSysx = kron(BlockI, BlockSysx);
        BlockSydx = kron(BlockSydx, BlockI);
        BlockSzsx = kron(BlockI, BlockSzsx);
        BlockSzdx = kron(BlockSzdx, BlockI);
        BlockI = kron(BlockI, BlockI);
        
        BlockH = trans(O)*H*O;
        BlockSxsx = trans(O)*BlockSxsx*O;
        BlockSxdx = trans(O)*BlockSxdx*O;
        BlockSysx = trans(O)*BlockSysx*O;
        BlockSydx = trans(O)*BlockSydx*O;
        BlockSzsx = trans(O)*BlockSzsx*O;
        BlockSzdx = trans(O)*BlockSzdx*O;
        BlockI = trans(O)*BlockI*O;
        
        cout << "Renormalization of sigmaZ..." << endl;
        for(int k=0; k<pow(2, i+1); k++)
        {
            sigmaZnewnew[k] = trans(O)*sigmaZnew[k]*O;
        }
        
        
        for(int k=0; k<pow(2,i+1)/2; k++)
        {
            CnewnewA[k] = trans(O)*CnewA[k]*O;
            DsliceA[k] = 0.5*gamma*( - kron(CnewnewA[k].st()*conj(CnewnewA[k]),BlockI) - kron(BlockI,trans(CnewnewA[k])*CnewnewA[k]) + 2.*kron(conj(CnewnewA[k]),CnewnewA[k]));
        }
        for(int k=0; k<pow(2,i+1)/2; k++)
        {
            CnewnewB[k] = trans(O)*CnewB[k]*O;
            DsliceB[k] = 0.5*gamma*( - kron(CnewnewB[k].st()*conj(CnewnewB[k]),BlockI) - kron(BlockI,trans(CnewnewB[k])*CnewnewB[k]) + 2.*kron(conj(CnewnewB[k]),CnewnewB[k]));
        }
        
        cx_mat GammaA = 0.*DsliceA[0];
        cx_mat GammaB = 0.*DsliceB[0];
        
        for(int k=0; k<pow(2,i+1)/2; k++)
        {
            GammaA = GammaA + DsliceA[k];
        }
        
        for(int k=0; k<pow(2,i+1)/2; k++)
        {
            GammaB = GammaB + DsliceB[k];
        }
        Da = GammaA;
        Db = GammaB;
        Ca = CnewnewA;
        Cb = CnewnewB;
        sigmaZ = sigmaZnewnew;
        
        cout << "Studying the liouvillian in the corner-space..." << endl;
        L_AUB = ii*kron(BlockH.st(), BlockI) - ii*kron(BlockI, BlockH) + Da + Db;
        LiouvSO L(L_AUB);
        dmAUB = L.GetDM();
        
        for(int k=0; k<pow(2, i+1); k++)
        {
            myfile << k << "\t" << real(trace(dmAUB*sigmaZ[k])) << endl;
            cout << "<sZ_totale[" << k << "]> = " << trace(dmAUB*sigmaZ[k]) << endl;
        }
        
        cout << "The " << i << " loop is finished and all seems to be alright" << endl;
    }
    cout << endl << endl;
    
    return 0;
}
