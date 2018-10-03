// Heisenberg chain
// boundary driven heisenberg chain

#include <iostream>
#include <armadillo>
#include <iomanip>
#include <complex>
#include <vector>
#include <fstream>
#include <cstdlib>
#include "init_matrices.hpp"
#include "HMatrix.hpp"
#include "utility.hpp"

using namespace std;
using namespace arma;
using namespace corner;

int main(int argc, char *argv[])
{
    
    double gamma = 1.;
    int Niter = atoi(argv[1]);
    int MaxM = atoi(argv[2]);
    float Jx = 1.; //atof(argv[3]);
    float Jy = 0.5; //atof(argv[4]);
    ofstream myfile("sz_HeisChain.txt");
    
    int M = MaxM;
    //for(int M=1; M<=MaxM; M++)
    //{
        cx_mat BlockI, I; //identity matrix
        cx_double ii; //imaginary unit
        cx_mat b; //jump operator
        cx_mat Sx, Sy, Sz, BlockH; //H's matrices
        init_matrices(Sx, Sy, Sz, I, b, ii, BlockH);
        cx_mat BlockSxsx = Sx;
        cx_mat BlockSxdx = Sx;
        cx_mat BlockSysx = Sy;
        cx_mat BlockSydx = Sy;
        cx_mat BlockSz = Sz;
        cx_mat H = BlockH;
        BlockI = I;
        
        cx_mat D;
        vector<cx_mat> C(pow(2,Niter));
        vector<cx_mat> Cnew(pow(2,Niter));
        vector<cx_mat> Cnewnew(pow(2,Niter));
        vector<cx_mat> Dslice(pow(2,Niter));
        
        C[0] = b;
        
        D = -0.5*gamma*kron(C[0].st()*conj(C[0]),BlockI) - 0.5*gamma*kron(BlockI,trans(C[0])*C[0]) + kron(conj(C[0]),C[0]);
        
        cx_mat L; //liouvillian superoperator
        
        for(int i=0; i<Niter; i++)
        {
            L = ii*kron(BlockH.st(), BlockI) - ii*kron(BlockI, BlockH) + D;
            //cout << "Ldim: " << L.n_rows << "x" << L.n_cols << endl;
            //cout << "Idim: " << BlockI.n_rows << "x" << BlockI.n_cols << endl;
            //cout << "Hdim: " << BlockH.n_rows << "x" << BlockH.n_cols << endl;
            
            cx_mat Evect;
            cx_vec Spectrum;
            eig_gen(Spectrum, Evect, L);
            
            cx_vec SortedSpectrum = sort(Spectrum);
            uvec indicesSpectrum = sort_index(Spectrum);
            cx_mat SortedEvect = Evect.cols(indicesSpectrum);
            
            cout << "Spectrum of L: " << endl << fixed << setprecision(10) << SortedSpectrum << endl;
            
            cx_vec rho_ss = SortedEvect.col(0);
            
            cx_mat dm = reshape(rho_ss, sqrt(L.n_rows), sqrt(L.n_rows));
            dm = dm/trace(dm);
            //cout << "Tr(dm) = " << trace(dm) << endl;// << "dm: " << dm << endl;
            bool equal = approx_equal(trans(dm), dm, "absdiff", 1E-10);
            //cout << "Is the dm hermitian? " << equal << endl;
            
            HMatrix DM(dm);
            cx_mat EvectDm = DM.GetONBasis();
            cx_vec SpectrumDm = DM.GetEigenvalues();
            
            //cout << "DM eigval: " << endl << SpectrumDm << endl;
            
            if(i == Niter-1)
            {
                cx_double expValue = trace(dm*BlockSz);
                myfile << M << "\t" << real(expValue) << "\n";
            }
            
            cx_vec SpectrumAUB = kron(SpectrumDm, SpectrumDm);
            cx_mat EvectAUB = kron(EvectDm, EvectDm);
            
            cx_vec SortedSpectrumAUB = sort(SpectrumAUB);
            uvec indices_aub = sort_index(SpectrumAUB);
            cx_mat SortedEvectAUB = EvectAUB.cols(indices_aub);
            
            int EvectCols = SortedEvectAUB.n_cols;
            vec choice;
            choice << EvectCols << M << endr;
            int Nkeep = min(choice); //valore minimo del vettore choice
            
            cx_mat O = SortedEvectAUB.tail_cols(Nkeep); //matrice pseudo-ortogonale
            
            //cout << "Odim: " << O.n_rows << "x" << O.n_cols << endl;
            
            for(int k=0; k<pow(2,i+1); k++)
            {
                if(k<pow(2,i))
                {
                    Cnew[k]=kron(BlockI, C[k]);
                }
                else
                {
                    Cnew[k] = kron(BlockI, C[k-pow(2,i)]);
                }
            }
            
            H = kron(BlockH, BlockI) + kron(BlockI, BlockH) - Jx*kron(BlockSxsx, BlockSxdx) - Jy*kron(BlockSysx, BlockSydx);
            
            BlockSxsx = kron(BlockI, BlockSxsx);
            BlockSxdx = kron(BlockSxdx, BlockI);
            BlockSysx = kron(BlockI, BlockSysx);
            BlockSydx = kron(BlockSydx, BlockI);
            BlockSz = kron(BlockI, BlockSz);
            BlockI = kron(BlockI, BlockI);
            //cout << "Idim: " << BlockI.n_rows << "x" << BlockI.n_cols << endl;
            
            BlockH = trans(O)*H*O;
            BlockSxsx = trans(O)*BlockSxsx*O;
            BlockSxdx = trans(O)*BlockSxdx*O;
            BlockSysx = trans(O)*BlockSysx*O;
            BlockSydx = trans(O)*BlockSydx*O;
            BlockSz = trans(O)*BlockSz*O;
            BlockI = trans(O)*BlockI*O;
            
            //cout << "Ddim: " << D.n_rows << "x" << D.n_cols << endl;
            //cout << "Cnewnewdim: " << Cnewnew[0].n_rows << "x" << Cnewnew[0].n_cols << endl;
            //cout << "Cnewdim: " << Cnew[0].n_rows << "x" << Cnew[0].n_cols << endl;
            
            for(int k=0; k<pow(2,i+1); k++)
            {
                Cnewnew[k] = trans(O)*Cnew[k]*O;
                //cout << "Cnewdim: " << Cnew[0].n_rows << "x" << Cnew[0].n_cols << endl;
                //cout << "Cnewnewdim: " << Cnewnew[0].n_rows << "x" << Cnewnew[0].n_cols << endl;
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
            //cout << "Ddim: " << D.n_rows << "x" << D.n_cols << endl;
        }
        //cout << endl << endl;
    //}
    
    return 0;
}
