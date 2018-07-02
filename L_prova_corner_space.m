%Heisenberg model, prova: L per quattro siti
clear all

gamma = 1.;
Delta = 0.;

BlockI = eye(2);
Sz = [1 0; 0 -1];
Sx = [0 1; 1 0];
Sy = [0 -1i; 1i 0];

BlockSz = Sz;
BlockSx = Sx;
BlockSy = Sy;

%jump operator
b = [0 0; 2 0];
BlockC = b;

H = - kron(Sx, kron(Sx, kron(BlockI, BlockI))) ...
    - kron(BlockI, kron(Sx, kron(Sx, BlockI))) ...
    - kron(BlockI, kron(BlockI, kron(Sx, Sx))) ...
    - 1/2*kron(Sy, kron(Sy, kron(BlockI, BlockI))) ...
    - 1/2*kron(BlockI, kron(Sy, kron(Sy, BlockI))) ...
    - 1/2*kron(BlockI, kron(BlockI, kron(Sy, Sy)));

C = kron(b, kron(BlockI, kron(BlockI, BlockI))) ...
    + kron(BlockI, kron(b, kron(BlockI, BlockI))) ...
    + kron(BlockI, kron(BlockI, kron(b, BlockI))) ...
    + kron(BlockI, kron(BlockI, kron(BlockI, b)));

BlockI = kron(BlockI, kron(BlockI, kron(BlockI, BlockI)));

%H' è la trasposta di H
first = 1i * H' - gamma * 1/2 * C.' * conj(C);
second = -1i * H - gamma * 1/2 * C' * C;
third = gamma * kron(conj(C), C);

L = kron(first, BlockI) + kron(BlockI, second) + third;

[Evect Spectrum] = eig(L);
Spectrum = diag(Spectrum);
SortedSpectrum = sort(Spectrum, 'ComparisonMethod', 'abs');

disp(SortedSpectrum);


