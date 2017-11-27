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
    int Niter = 1;
    double gamma = 0.5;
    double Delta = 0.;
    int n = 4; //number of sites
    
    mat I = eye<mat>(2,2);
    
    //definition of the imaginary unit
    complex <double> ii (0, 1.);
    
    //annihilation operator b's definition (fermionic case) - matrice 2x2
    mat b;
    b << 0 << 1 << endr << 0 << 0 << endr;
    
    
    //jumps operator's C definition - matrice 2x2
    mat C = gamma * b;
    
    
    //definition of the hamiltonian
    mat Sz, Sp, Sm;
    Sz << 0.5 << 0 << endr << 0 << -0.5 << endr;
    Sp << 0 << 0 << endr << 1 << 0 << endr;
    Sm << 0 << 1 << endr << 0 << 0 << endr;

    mat BlockSz = Sz;   //2x2
    mat BlockSp = Sp;   //2x2
    mat BlockSm = Sm;   //2x2
    mat BlockI = I;     //2x2
    mat BlockH = zeros<mat>(2,2);   //2x2
    mat h;
    double Energy = 0.;
    double LastEnergy;
    
    
    //inizializzazione dell'operatore di Lindblad
    cx_mat first = -ii * BlockH - gamma * trans(C) * C;  //primo termine di L
    
    cx_mat second = ii * trans(BlockH) - gamma * trans(C) * C;  //secondo termine di L
    
    mat third = gamma * kron(C, trans(C));  //terzo termine di L - 4x4
    
    cx_mat L = kron(first, I) + kron(I, second) + third;  //L - matrice 4x4
    L = 0.5 * (L + trans(L));
    
    LastEnergy = Energy;
    
    
    //diagonalizzazione di L
    vec Spectrum;
    cx_mat Evect;
    eig_sym(Spectrum, Evect, L);
    int i = 0;
    
    double SystSize, Ener2, EnergyPerBond;
    SystSize = pow(2,i);
    
    Energy = Spectrum[0];
    Ener2 = Energy/SystSize;
    EnergyPerBond = (Energy-LastEnergy)/(SystSize/2);
    
     
    //troncamento
    int Nkeep, SpectrumRows;
    cx_mat Omatr;  //matrice degli autovettori di L
     
    SpectrumRows = Spectrum.n_rows;  //#righe della matrice Spectrum (autovalori di L)
     
    Nkeep = min(SpectrumRows, m);  //prendo il minimo tra #righe di Spectrum e m (scelto da me)
     
    Omatr = Evect.cols(0, Nkeep-1);  //matrice formata dagli autovettori di L (presi fino a Nkeep)
     
    //cout << fixed << setprecision(6) << SystSize << "\t" << Energy << "\t" << EnergyPerBond << "\t" << Ener2 << endl;
    cout << Omatr << endl;
    
    mat O = real(Omatr);
    
    cout << O << endl;
    L = trans(Omatr)*L*Omatr;
    
    cout << "Dim(L): " << L.n_rows << " x " << L.n_cols << endl;
    
    //le seguenti diventano matrici 4x4
    BlockH = kron(BlockI, BlockH);
    //cout << "BlockH è una matrice " << BlockH.n_rows << " x " << BlockH.n_cols << endl;
    C = kron(BlockI, C);
    BlockSz = kron(BlockI, BlockSz);
    //cout << "BlockSz è una matrice " << BlockSz.n_rows << " x " << BlockSz.n_cols << endl;
    BlockSp = kron(BlockI, BlockSp);
    //cout << "BlockSp è una matrice " << BlockSp.n_rows << " x " << BlockSp.n_cols << endl;
    BlockSm = kron(BlockI, BlockSm);
    //cout << "BlockSm è una matrice " << BlockSm.n_rows << " x " << BlockSm.n_cols << endl;
    BlockI = kron(BlockI, BlockI);
    //cout << "BlockI è una matrice " << BlockI.n_rows << " x " << BlockI.n_cols << endl;
    

    //cambiamento di base
    BlockH = trans(O)*BlockH*O;
    C = trans(O)*C*O;
    BlockSz = trans(O)*BlockSz*O;
    BlockSp = trans(O)*BlockSp*O;
    BlockSm = trans(O)*BlockSm*O;
    BlockI = trans(O)*BlockI*O;
    
    //cout << L << endl;
    
    //for (int i=1; i<Niter+1; i++)
    //{
        int SystDim;
        SystDim = pow(2, 2*i);
        
        h = kron(BlockH, BlockI) + kron(BlockI, BlockH) - Delta * kron(BlockSz, BlockSz) + 0.5 * (kron(BlockSp, BlockSm) + kron(BlockSm, BlockSp));     //h = 4x4 matrix
        h = 0.5 * (h + trans(h));
    
    cout << "Dim(h): " << h.n_rows << " x " << h.n_cols << endl;
    
      /*
        //definition of the Lindblad operator
        
        cx_mat first = -ii * BlockH - gamma * trans(C) * C;
        
        cx_mat second = ii * trans(BlockH) - gamma * trans(C) * C;
        
        mat third = gamma * kron(C, trans(C));
        
        cx_mat L = kron(first, I) + kron(I, second) + third;
        L = 0.5 * (L + trans(L));
        
        LastEnergy = Energy;
        
        
        //diagonalization
        vec Spectrum;
        cx_mat Evect;
        eig_sym(Spectrum, Evect, L);
        
        cout << Evect << endl;
        
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
        
        L = trans(Omatr)*L*Omatr;
       
        BlockH = kron(BlockI, BlockH);
        BlockSz = kron(BlockI, BlockSz);
        BlockSp = kron(BlockI, BlockSp);
        BlockSm = kron(BlockI, BlockSm);
        BlockI = kron(BlockI, BlockI);
        
        
        BlockH = trans(Omatr)*h*Omatr;
        BlockSz = trans(Omatr)*BlockSz*Omatr;
        BlockSp = trans(Omatr)*BlockSp*Omatr;
        BlockSm = trans(Omatr)*BlockSm*Omatr;
        BlockI = trans(Omatr)*BlockI*Omatr;
    }
         */
         
    
    return 0;
}
