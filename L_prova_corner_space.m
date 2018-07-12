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
b = [0 0; 1 0];
BlockC = b;

H = - kron(Sx, kron(Sx, kron(BlockI, BlockI))) ...
    - kron(BlockI, kron(Sx, kron(Sx, BlockI))) ...
    - kron(BlockI, kron(BlockI, kron(Sx, Sx))) ...
    - 0.5*(kron(Sy, kron(Sy, kron(BlockI, BlockI))) ...
    +kron(BlockI, kron(Sy, kron(Sy, BlockI))) ...
    +kron(BlockI, kron(BlockI, kron(Sy, Sy))));

C1 = kron(b, kron(BlockI, kron(BlockI, BlockI)));
C2 = kron(BlockI, kron(b, kron(BlockI, BlockI)));
C3 = kron(BlockI, kron(BlockI, kron(b, BlockI)));
C4 = kron(BlockI, kron(BlockI, kron(BlockI, b)));

I = kron(BlockI, kron(BlockI, kron(BlockI, BlockI)));

%H' è la trasposta di H
first = 1i * H.' - gamma * 1/2 * C1.' * conj(C1) - gamma * 1/2 * C2.' * conj(C2) ...
    - gamma * 1/2 * C3.' * conj(C3) - gamma * 1/2 * C4.' * conj(C4);
second = -1i * H - gamma * 1/2 * (C1') * C1 - gamma * 1/2 * (C2') * C2 ...
    - gamma * 1/2 * (C3') * C3 - gamma * 1/2 * (C4') * C4;
third = gamma * kron(conj(C1), C1) + gamma * kron(conj(C2), C2) + ...
    + gamma * kron(conj(C3), C3) + gamma * kron(conj(C4), C4);

L = kron(first, I) + kron(I, second) + gamma * third;

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

%dmT = dm';
%equal = isequal(dmT, dm);

%if(equal)
%    disp("p è hermitiana.");
%end

if(abs(dm-dm') < exp(-14))
   disp("p è hermitiana");  %ordine di E-15
end


disp(trace(dm));

dens_matr = dm/trace(dm)
disp(trace(dens_matr));

Sz = kron(Sz, kron(BlockI, kron(BlockI, BlockI)));

expValueSz = trace(dens_matr * Sz)







