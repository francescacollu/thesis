%corner-space renormalization method
%Il metodo di ordinamento secondo gli indices è testato e funzionante.

gamma = 1.;
Delta = 0.;
Niter = 2;
M = 17;

%for M = 4:+4:44
    
    BlockI = eye(2);
    Sz = [1 0; 0 -1];
    Sx = [0 1; 1 0];
    Sy = [0 -1i; 1i 0];

    BlockSz = Sz;
    BlockSxsx = Sx;
    BlockSxdx = Sx;
    BlockSysx = Sy;
    BlockSydx = Sy;

    %jump operator
    b = [0 0; 1 0];
    C = b;

    BlockH = zeros(2);
    H = BlockH;

    
    for k = 1:Niter

        first = 1i * (BlockH.') - 0.5 * (C.' * conj(C)); %C.' = trasposto di C; H' = trasposta compl.coniugata di H
        first = kron(first, BlockI);
        second = -1i * (BlockH) - 0.5 * (C' * C);
        second = kron(BlockI, second);
        third = kron(conj(C), C);

        L = first + second + third;

        [Evect, Spectrum] = eig(L);
        [SortedSpectrum, Indices] = sort(diag(Spectrum), 'ComparisonMethod', 'abs');
        EvectSorted = Evect(:, Indices); %ordinamento di Evect secondo gli indici di ordinamento di Spectrum - testato, funzionante

        disp(SortedSpectrum);

        rhoss = EvectSorted(:, 1);
        szDimL = size(L, 1);
        sz = sqrt(szDimL);
        rho = reshape(rhoss, [sz, sz]);
        dm = rho;

        if(k == Niter) 
            display("M : "+ M + "   Tr(p*SigmaZ) = "+ trace(dm*BlockSz));
        end

        % Merging
        [eigenvecDM, spectrumDM] = eig(dm);
        eigenvecDmAUB = kron(eigenvecDM, eigenvecDM);
        spectrumDmAUB = kron(spectrumDM, spectrumDM);

        [spectrumDmAUBsorted, Index] = sort(diag(spectrumDmAUB), 'ComparisonMethod', 'abs');
        eigenvecDmAUBsorted = eigenvecDmAUB(:, Index);

        Z = [size(eigenvecDmAUBsorted, 2); M];
        Nkeep = min(Z); %valore minimo del vettore Z
        Omatr = eigenvecDmAUBsorted(:, end-Nkeep+1:end);

        H = kron(BlockH, BlockI) + kron(BlockI, BlockH) - kron(BlockSxsx, BlockSxdx) ...
        - 0.5 * (kron(BlockSysx, BlockSydx));

        C = kron(C, BlockI) + kron(BlockI, C);
        
        BlockSxsx = kron(BlockI, BlockSxsx);
        BlockSxdx = kron(BlockSxdx, BlockI);
        BlockSysx = kron(BlockI, BlockSysx);
        BlockSydx = kron(BlockSydx, BlockI);
        BlockSz = kron(BlockI, BlockSz);
        BlockI = kron(BlockI, BlockI);

        BlockH = (Omatr)'*H*Omatr;
        C = (Omatr)'*C*Omatr;
        BlockSxsx = (Omatr)'*BlockSxsx*Omatr;
        BlockSxdx = (Omatr)'*BlockSxdx*Omatr;
        BlockSysx = (Omatr)'*BlockSysx*Omatr;
        BlockSydx = (Omatr)'*BlockSydx*Omatr;
        BlockSz = (Omatr)'*BlockSz*Omatr;
        BlockI = (Omatr)'*BlockI*Omatr;


    end

