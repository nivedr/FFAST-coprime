
function [DFT] = FFAST (FFTarray, N1, N2, d, sampleArray)

% N1, N2 -> dimensions of input array
% d -> number of delays
% sampleArray -> array with the number of times to subsample

    array = ifft2(FFTarray);

    dim = N1*N2;
    K = zeros(N1,N2);
    N = zeros(N1,N2);
    
    for n = (0 : N1-1);
        for m = (0 : N2-1);
            K(n+1,m+1) = mod(N2*n + N1*m, dim);
        end
    end
    K = reshape(K,1,N1*N2);
    
    for k = (0 : N1-1)
        for l = (0 : N2-1)
            N(k+1,l+1) = mod(k*invMod(N2, N1)*N2 + l*invMod(N1, N2)*N1, dim);
        end
    end
    N = reshape(N,1,N1*N2);
    
    % q = [0, 16, 12, 8, 4; 5, 1, 17, 13, 9; 10, 6, 2, 18, 14; 15, 11, 7, 3, 19];
    
    p = (1:N1*N2);
    CRTArray(p) = array(N(p)+1);
    
    % cellArray = cell(length(sampleArray));
    
    for i = (1 : length(sampleArray))
        SubsampledArray = downsample(CRTArray(p), sampleArray(i));
        if d > 1
            for j = (1 : d-1)
                sample = downsample(CRTArray(p), sampleArray(i), j);
                SubsampledArray = [ SubsampledArray ; sample ];
            end
        end
        c = SubsampledArray.';
        SubsampledArrayFFT = sampleArray(i)*fft(c).';
        subsampledCell{i} = SubsampledArray;
        FFTsubsampledCell{i} = SubsampledArrayFFT;
    end
    % cellArray now holds in each cell a 2d array of subsampled and delayed
    % inputs
    
    DFT = PeelingSolver(FFTsubsampledCell, d, sampleArray, dim);
    %DFT = reshape(DFT, N2, N1)';
    finalDFT = zeros(N1, N2);

    [G, U, V] = gcd(N1, N2);    % G = 1, U & V solve Diophatine Equations
    for i = (0:dim-1)
        n1 = mod(V*i, N1);
        n2 = mod(U*i, N2);
        finalDFT(n1+1, n2+1) = DFT(i+1);
    end
    a = 1;
end
