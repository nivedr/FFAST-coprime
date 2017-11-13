
function [DFT] = FFAST (FFTarray, N1, N2, d, N1sampleArray, N2sampleArray)

% N1, N2 -> dimensions of input array
% d -> number of delays
% N1sampleArray, N2sampleArray :: same length

    array = ifft2(FFTarray);

    dim = N1*N2;
    K = zeros(N1,N2);
    N = zeros(N1,N2);    
    p = (1:N1*N2);
    
    % cellArray = cell(length(sampleArray));
    for i = (1 : length(N1sampleArray))
        sample00 = downsample(downsample(array(:,:), N1sampleArray(i))', N2sampleArray(i))';
        if d > 1
            for j = (1 : d-1)
                sample10 = downsample(downsample(array(:,:), N1sampleArray(i), 1)', N2sampleArray(i), 0)';
                sample01 = downsample(downsample(array(:,:), N1sampleArray(i), 0)', N2sampleArray(i), 1)';

                SubsampledArray = [ sample00 ; sample10 ; sample01 ];
            end
        end
        SubsampledArrayFFT = N2sampleArray(i)*N1sampleArray(i)*[ fft2(sample00) ; fft2(sample10) ; fft2(sample01) ];
        subsampledCell{i} = SubsampledArray;
        FFTsubsampledCell{i} = SubsampledArrayFFT;
        N1stride = N1/N1sampleArray(i);
        N2stride = N2/N2sampleArray(i);

        for j = (1:N1stride)
            for k = (1:N2stride) % indexing with stream no. also
                Yijk = [ SubsampledArrayFFT(j,k) ; SubsampledArrayFFT(j+N1stride,k) ; SubsampledArrayFFT(j+2*N1stride,k) ];
                Y{j,k} = Yijk;
            end
        end
        YStream{i} = Y;
        Y = [];
    end
    % cellArray now holds in each cell a 2d array of subsampled and delayed
    % inputs
    
    DFT = PeelingSolver(YStream, N1, N2, N1sampleArray, N2sampleArray);
    a = 1;
end
