//NRG
// H = -\sum [Sx Sx + Sy Sy + Delta Sz Sz]

#include <iostream>
#include <armadillo>
#include <iomanip>
#include <chrono>
#include <ctime>

using namespace std;
using namespace arma;

int main()
{
    auto start = chrono::system_clock::now();
    
    // set z anisotropy = 0
    double Delta = 0.;
    int m = 1;
    int Niter = 10;
    
    mat I = eye<mat>(2,2);
    
    mat Sz, Sp, Sm;
    Sz << 0.5 << 0 << endr << 0 << -0.5 << endr;

    Sp << 0 << 0 << endr << 1 << 0 << endr;

    Sm << 0 << 1 << endr << 0 << 0 << endr;
    
    //blocks initialization
    mat BlockSz = Sz;
    mat BlockSp = Sp;
    mat BlockSm = Sm;
    mat BlockI = I;
    mat BlockH = zeros<mat>(2,2);
    double Energy = 0.;
    double LastEnergy;
    
    //NRG
    
    for (int i=1; i<Niter+1; i++)
    {
        double SystSize, Ener2, EnergyPerBond;
        SystSize = pow(2,i);
        mat H_super;
        
        H_super = kron(BlockH, BlockI) + kron(BlockI, BlockH) - Delta * kron(BlockSz, BlockSz) + 0.5 * (kron(BlockSp, BlockSm) + kron(BlockSm, BlockSp));
        H_super = 0.5 * (H_super + trans(H_super));
        
        LastEnergy = Energy;
        
        //diagonalization
        vec Spectrum;
        mat Evect;
        eig_sym(Spectrum, Evect, H_super);
        
        Energy = Spectrum[0];
        Ener2 = Energy/SystSize;
        EnergyPerBond = (Energy-LastEnergy)/(SystSize/2); //???
        
        
        //truncation operator
        int Nkeep, SpectrumRows;
        mat Omatr;
        
        //cout << size(Spectrum) << endl;
        SpectrumRows = Spectrum.n_rows;
        
        Nkeep = min(SpectrumRows, m);
        //cout << Nkeep << endl;
        Omatr = Evect.cols(0, Nkeep-1);
        
        cout << fixed << setprecision(9) << SystSize << "\t" << Energy << "\t" << EnergyPerBond << "\t" << Ener2 << endl;
        
        BlockSz = kron(BlockI, BlockSz);
        BlockSp = kron(BlockI, BlockSp);
        BlockSm = kron(BlockI, BlockSm);
        BlockI = kron(BlockI, BlockI);
        
        BlockH = trans(Omatr)*H_super*Omatr;
        BlockSz = trans(Omatr)*BlockSz*Omatr;
        BlockSp = trans(Omatr)*BlockSp*Omatr;
        BlockSm = trans(Omatr)*BlockSm*Omatr;
        BlockI = trans(Omatr)*BlockI*Omatr;
        
    }
    
    auto end = chrono::system_clock::now();
    
    chrono::duration<double> elapsed_seconds = end-start;
    
    time_t end_time = chrono::system_clock::to_time_t(end);
    
    cout << "Finished computation at " << ctime(&end_time) << endl;
    cout << "Elapsed time: " << elapsed_seconds.count() << "s" << endl;
    
    return 0;
}
