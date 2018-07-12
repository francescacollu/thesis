%Heisenberg model, prova: L per due siti
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
b = [0 0; 1 0];
C1 = kron(b, BlockI);
C2 = kron(BlockI, b);

I = kron(BlockI, BlockI);

H = - kron(Sx, Sx) - 0.5*kron(Sy, Sy);

%H' è la trasposta di H
first = 1i * H' - gamma * 1/2 * C1.' * conj(C1) - gamma * 1/2 * C2.' * conj(C2);
second = -1i * H - gamma * 1/2 * (C1') * C1 - gamma * 1/2 * (C2') * C2;
third = gamma * kron(conj(C1), C1) + gamma * kron(conj(C2), C2);

L = kron(first, I) + kron(I, second) + third;

[Evect, Spectrum] = eig(L);
[SortedSpectrum, Indices] = sort(diag(Spectrum), 'ComparisonMethod', 'abs');
EvectSorted = Evect(:, Indices);

disp(SortedSpectrum);

rhoss = EvectSorted(:, 1);
szDimL = size(L, 1);
sz = sqrt(szDimL);
rho = reshape(rhoss, [sz, sz]);
dm = rho;
[Edm, SpectrumDM] = eig(dm);
display(diag(SpectrumDM));

display("Somma(lambda_i) = " + sum(diag(SpectrumDM)));

dmT = dm';
equal = isequal(dmT, dm);

if(equal)
    disp("p è hermitiana.");
end

disp(trace(dm));

dens_matr = dm/trace(dm)
disp(trace(dens_matr));

Sz = kron(Sz, BlockI);

expValueSz = trace(dens_matr * Sz)




