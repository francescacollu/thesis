%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%               NRG for the  1D  XXZ  Heisenberg Model
%              H = - \sum [ Sx Sx + Sy Sy + Delta Sz Sz ]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all

tic
Delta = 0.;   % ZZ anisotropy
m=10;          % Number of states kept m
NIter = 10;    % Number of iterations.  Final lattice size is 2*Niter + 2


%               Intialize local operators
I= eye(2);
Sz = [1/2 0 ; 0 -1/2];
Sp = [0 0 ; 1 0];
Sm = [0 1 ; 0 0];

% Initial blocks  (we assume reflection symmetry)
BlockSz = Sz;
BlockSp = Sp;
BlockSm = Sm;
BlockI  = I;
BlockH  = zeros(2);
Energy = 0.;

%%  Numerical Renormalization Group (NRG)

for l = 1:NIter
    
    SystSize = 2^l;
    
%                   HAMILTONIAN MATRIX forsuperblock
    H_super = kron(BlockH, BlockI) + kron(BlockI, BlockH) ...
        - Delta * kron(BlockSz, BlockSz) ...
        + 0.5 * ( kron(BlockSp, BlockSm) + kron(BlockSm, BlockSp) );

    H_super = 0.5 * (H_super + H_super');  % ensure H is symmetric

    LastEnergy = Energy;
    
%                   Diagonalizing the Hamiltonian
    [Evect Spectrum] = eig(H_super);
    [Spectrum Index] = sort(diag(Spectrum), 'ascend');  % sort ascending
    Evect = Evect(:,Index);
    
    Energy = Spectrum(1);
    Ener2  = Energy / SystSize;
    EnergyPerBond = (Energy - LastEnergy) / (SystSize/2);
    
%                   Construct the truncation operator
    NKeep = min(size(Spectrum, 1), m);
    Omatr = Evect(:, 1:NKeep);

    fprintf('%d\t%f\t%f\t%f\n',SystSize, Energy, EnergyPerBond, Ener2);

    BlockSz = kron(BlockI,BlockSz);
    BlockSp = kron(BlockI,BlockSp);
    BlockSm = kron(BlockI,BlockSm);
    BlockI = kron(BlockI,BlockI);
    
%  Transform the block operators into the truncated basis
    BlockH  = Omatr' * H_super * Omatr;
    BlockSz = Omatr' * BlockSz * Omatr;
    BlockSp = Omatr' * BlockSp * Omatr;
    BlockSm = Omatr' * BlockSm * Omatr;
    BlockI =  Omatr' * BlockI  * Omatr;

end

toc
