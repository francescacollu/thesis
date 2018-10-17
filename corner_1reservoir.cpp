// corner-space-renormalization-method
// codice strutturato in classi

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
    int M = 30;
    ofstream myfile("expValueSz.txt");

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
    cx_mat BlockSz = Sz;
    cx_mat H = BlockH;
    BlockI = I;

    vector<cx_mat> sigmaZ(pow(2, Niter));
    vector<cx_mat> sigmaZnew(pow(2, Niter));
    vector<cx_mat> sigmaZnewnew(pow(2, Niter));

    sigmaZ[0] = Sz;
    
    cx_mat D;
    vector<cx_mat> C(pow(2,Niter));
    vector<cx_mat> Cnew(pow(2,Niter));
    vector<cx_mat> Cnewnew(pow(2,Niter));
    vector<cx_mat> Dslice(pow(2,Niter));
    
    C[0] = b;
    
    D = -0.5*gamma*kron(C[0].st()*conj(C[0]),BlockI) - 0.5*gamma*kron(BlockI,trans(C[0])*C[0]) + kron(conj(C[0]),C[0]);
    //Db = -0.5*gamma*kron(BlockI.st()*conj(BlockI),BlockI) - 0.5*gamma*kron(BlockI,trans(BlockI)*BlockI) + kron(conj(BlockI),BlockI);
    
    cx_mat L_A, L_B; //liouvillian superoperator
    
    for(int i=0; i<Niter; i++)
    {
        ////////// First block //////////
        L_A = ii*kron(BlockH.st(), BlockI) - ii*kron(BlockI, BlockH) + D;
        cout << "L_Adim: " << L_A.n_rows << "x" << L_A.n_cols << endl;
        cout << "***DM A***" << endl;
        LiouvSO LA(L_A);
        cx_mat dmA = LA.GetDM();
        HMatrix DM_A(dmA);
        cx_mat EvectDmA = DM_A.GetONBasis();
        cx_vec SpectrumDmA = DM_A.GetEigenvalues();
        
        ////////// Second block //////////
        L_B = ii*kron(BlockH.st(), BlockI) - ii*kron(BlockI, BlockH);
        cout << "L_Bdim: " << L_B.n_rows << "x" << L_B.n_cols << endl;
        cout << "***DM B***" << endl;
        LiouvSO LB(L_B);
        cx_mat dmB = LB.GetDM();
        HMatrix DM_B(dmB);
        cx_mat EvectDmB = DM_B.GetONBasis();
        cx_vec SpectrumDmB = DM_B.GetEigenvalues();
        
        //cout << "Idim: " << BlockI.n_rows << "x" << BlockI.n_cols << endl;
        //cout << "Hdim: " << BlockH.n_rows << "x" << BlockH.n_cols << endl;
        
        bool test = approx_equal(dmB*EvectDmB.col(dmB.n_rows - 1), EvectDmB.col(dmB.n_rows - 1)*SpectrumDmB(dmB.n_rows - 1), "absdiff", 1E-10);
        
        cout << "Is there a match? " << test << endl;
        
        //cout << "DM eigval: " << endl << SpectrumDm << endl;
        
        if(i == Niter-1)
        {
            for(int k=0; k<pow(2, i+1)/2; k++)
            {
                myfile << k << "\t" << real(trace(dmA*sigmaZ[k])) << endl;
                cout << "<sigmaZ>[" << k << "] = " << trace(dmA*sigmaZ[k]) << endl;
            }
        }
        
        ////////// MERGING the two blocks /////////
        cx_mat dmAUB = kron(dmA, dmB);
        HMatrix DM_AUB(dmAUB);
        cx_mat EvectDmAUB = DM_AUB.GetONBasis();
        cx_vec SpectrumDmAUB = DM_AUB.GetEigenvalues();
        //cout << "dmRows = " << dmAUB.n_rows << "x" << "dmCols = " << dmAUB.n_cols << endl;
        
        int EvectCols = EvectDmAUB.n_cols;
        vec choice;
        choice << EvectCols << M << endr;
        int Nkeep = min(choice); //valore minimo del vettore choice
        
        cx_mat O = EvectDmAUB.tail_cols(Nkeep); //matrice pseudo-ortogonale
        //cout << "Odim: " << O.n_rows << "x" << O.n_cols << endl;
        
        /*cout << "Idim: " << BlockI.n_rows << "x" << BlockI.n_cols << endl;
        cout << "Sxdim: " << BlockSxdx.n_rows << "x" << BlockSxdx.n_cols << endl;
        cout <<  BlockSxsx.n_rows << "x" << BlockSxsx.n_cols << endl;
        cout << "Sydim: " << BlockSydx.n_rows << "x" << BlockSydx.n_cols << endl;
        cout << BlockSysx.n_rows << "x" << BlockSysx.n_cols << endl;
        cout << "Szdim: " << BlockSzdx.n_rows << "x" << BlockSzdx.n_cols << endl;
        cout << BlockSzsx.n_rows << "x" << BlockSzsx.n_cols << endl;
        cout << "Hdim: " << BlockH.n_rows << "x" << BlockH.n_cols << endl;*/
        
        
        ////////// Merging the Hamiltonians //////////
        H = kron(BlockH, BlockI) + kron(BlockI, BlockH) - Jx*kron(BlockSxsx, BlockSxdx) - Jy*kron(BlockSysx, BlockSydx) - Jz*kron(BlockSzsx, BlockSzdx);
        
        // only the first block has (one) dissipator (in the 0 site)
        for(int k=0; k<pow(2, i+1); k++) // distribution of the dissipator(s) over the first block-for the next loop
        {
            if(k<pow(2,i))
            {
                Cnew[k]=kron(C[k], BlockI);
            }
            else
            {
                Cnew[k] = kron(BlockI, BlockI);//C[k-pow(2,i)]);
            }
        }
        
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
        BlockSz = kron(BlockSz, BlockI);
        BlockI = kron(BlockI, BlockI);
        //cout << "Idim: " << BlockI.n_rows << "x" << BlockI.n_cols << endl;
        
        //cout << "BlockIdim: " << BlockI.n_rows << "x" << BlockI.n_cols << endl;
        //cout << "SigmaZdim: " << sigmaZnew[0].n_rows << "x" << sigmaZnew[0].n_cols << endl;
        //cout << "DM_AUBdim: " << dmAUB.n_rows << "x" << dmAUB.n_cols << endl;
        
        /*if(i == Niter-1)
        {
            for(int k=0; k<pow(2, i+1); k++)
            {
                myfile << k << "\t" << real(trace(dmAUB*sigmaZnew[k])) << endl;
                cout << "<sigmaZ>[" << k << "] = " << trace(dmAUB*sigmaZnew[k]) << endl;
            }
        }*/
        
        BlockH = trans(O)*H*O;
        BlockSxsx = trans(O)*BlockSxsx*O;
        BlockSxdx = trans(O)*BlockSxdx*O;
        BlockSysx = trans(O)*BlockSysx*O;
        BlockSydx = trans(O)*BlockSydx*O;
        BlockSzsx = trans(O)*BlockSzsx*O;
        BlockSzdx = trans(O)*BlockSzdx*O;
        BlockSz = trans(O)*BlockSz*O;
        BlockI = trans(O)*BlockI*O;
        
        for(int k=0; k<pow(2, i+1); k++)
        {
            sigmaZnewnew[k] = trans(O)*sigmaZnew[k]*O;
        }
        
        //cout << "Ddim: " << D.n_rows << "x" << D.n_cols << endl;
        //cout << "Cnewnewdim: " << Cnewnew[0].n_rows << "x" << Cnewnew[0].n_cols << endl;
        //cout << "Cnewdim: " << Cnew[0].n_rows << "x" << Cnew[0].n_cols << endl;
        for(int k=0; k<pow(2,i+1); k++)
        {
            Cnewnew[k] = trans(O)*Cnew[k]*O;
            //cout << "Cnewdim: " << Cnew[k].n_rows << "x" << Cnew[k].n_cols << endl;
            //cout << "Cnewnewdim: " << Cnewnew[k].n_rows << "x" << Cnewnew[k].n_cols << endl;
            //cout << "BlockIdim: " << BlockI.n_rows << "x" << BlockI.n_cols << endl;
            Dslice[k] = 0.5*gamma*( - kron(Cnewnew[k].st()*conj(Cnewnew[k]),BlockI) - kron(BlockI,trans(Cnewnew[k])*Cnewnew[k]) + 2.*kron(conj(Cnewnew[k]),Cnewnew[k]));
            //cout << "Dslicedim: " << Dslice[0].n_rows << "x" << Dslice[0].n_cols << endl;
        }
        
        cx_mat Gamma = 0.*Dslice[0];
        for(int k=0; k<pow(2,i+1); k++)
        {
            Gamma = Gamma + Dslice[k];
            //cout << "Dslicedim: " << Dslice[k].n_rows << "x" << Dslice[k].n_cols << endl;
        }
        D = Gamma;
        C = Cnewnew;
        sigmaZ = sigmaZnewnew;
        //cout << "Ddim: " << D.n_rows << "x" << D.n_cols << endl;
        
        cout << "The " << i << " loop is finished and all seems to be alright" << endl;
    }
    cout << endl << endl;
    
    return 0;
}
