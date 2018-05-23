%corner-space renormalization method

clear all

gamma = 1.;
Delta = 0.;
Niter = 3;
M = 17;

BlockI = eye(2);
Sz = [1 0; 0 -1];
Sx = [0 1; 1 0];
Sy = [0 -1*1i; 1*1i 0];

BlockSz = Sz;
BlockSxsx = Sx;
BlockSxdx = Sx;
BlockSysx = Sy;
BlockSydx = Sy;

%jump operator
b = [0 0; 2 0];
C = b;

BlockH = zeros(2);
H = BlockH;

dm = zeros(2);

for l = 1:Niter
    
    first = 1i * (H') - gamma * 1/2 * (C' * conj(C));
    second = -1i * (H) - gamma * 1/2 * (C' * C);
    third = gamma * kron(conj(C), C);
    
    L = kron(first, BlockI) + kron(BlockI, second) + third;
    szDimL = size(L, 1);
    
    [Evect, Spectrum] = eig(L);
    [Spectrum, Indices] = sort(diag(Spectrum), 'ComparisonMethod', 'abs');
    
    disp(Spectrum);
    %fprintf('%f+%fi\n\n', real(Spectrum), imag(Spectrum));
    %fprintf('%f\n', size(Spectrum, 1));
    
    %{for k = 1:size(Spectrum, 1)
        %{if (abs(Spectrum(k)) < 1.E-10)
            %{rhoss = Evect(:, k);
            %szDimL = size(L, 1);
           % sz = sqrt(szDimL);
            %rho = reshape(rhoss, [sz, sz]);
                %if trace(rho) < 1.1
               % dm = rho;
                %fprintf (formatSpec, dm);
                %end
%end
    
    
    rhoss = Evect(:, 1);
    szDimL = size(L, 1);
    sz = sqrt(szDimL);
    rho = reshape(rhoss, [sz, sz]);
    dm = rho;
    [EvectRo, EvalRo] = eig(dm);
    lambdaVec = kron(EvalRo, EvalRo);
    phiVec = kron(EvectRo, EvectRo);
    
    %TrSz = trace(dm*BlockSz);
    %fprintf('%f\n', TrSz);
            
    %n=0;
            
    %for i = 1:size(EvectRo, 2)
    %    for j = 1:size(EvectRo, 2)
    %        eveci = Evect(:, i);
    %        evecj = Evect(:, j);
            
    %        phiphi = kron(eveci, evecj);
            
     %       phiVec(:, n) = phiphi;
     %       n = n+1;
     %   end
    %end
    
    
    [lambdaVecSorted, Index] = sort(diag(lambdaVec), 'ComparisonMethod', 'abs');
    phiVecSorted = phiVec(:, Index);
    
    Z = [size(phiVecSorted, 2); M];
    Nkeep = min(Z);
    Omatr = phiVecSorted(:, end-Nkeep+1:end);
    
    BlockH = H;
    
    H = kron(BlockH, BlockI) + ...
    kron(BlockI, BlockH) - kron(BlockSxsx, BlockSxdx) ...
    - 0.5 * kron(BlockSysx, BlockSydx);
    
    %disp(H);
    %disp(Omatr);
    
    BlockSz = kron(BlockI, BlockSz);
    BlockSxsx = kron(BlockI, BlockSxsx);
    BlockSxdx = kron(BlockSxdx, BlockI);
    BlockSysx = kron(BlockI, BlockSysx);
    BlockSydx = kron(BlockSydx, BlockI);
    C = kron(C, BlockI) + kron(BlockI, C);
    BlockI = kron(BlockI, BlockI);
    
    
    H = Omatr' * H * Omatr;
    C = Omatr' * C * Omatr;
    BlockSz = Omatr' * BlockSz * Omatr;
    BlockSxsx = Omatr' * BlockSxsx * Omatr;
    BlockSxdx = Omatr' * BlockSxdx * Omatr;
    BlockSysx = Omatr' * BlockSysx * Omatr;
    BlockSydx = Omatr' * BlockSydx * Omatr;
    BlockI = Omatr' * BlockI * Omatr;
    

end
    

