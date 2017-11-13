% len -> number of streams
% d -> number of sub-streams

function [DFT] = PeelingSolver(FFTsubsampledCell, d, sampleArray, dftSize)
    
    DFT = zeros(1, dftSize);
    len = length(sampleArray);
    
    empty = true;
    for i = (1 : len)
        FFTsubsampledCell{i} = round(FFTsubsampledCell{i},8);
        if any(any(FFTsubsampledCell{i})) ~= 0
            empty = false;
        end
    end
    
    if empty
        return;
    end
    
    for i = (1 : len)
        CurrentStream = FFTsubsampledCell{i};
        streamSize = sampleArray(i);
        for j = (1 : dftSize/streamSize)
            if CurrentStream(1,j) ~= 0
                Position = round((dftSize/(2*pi))*angle(CurrentStream(2,j)*conj(CurrentStream(1,j))), 4);
            else
                % The node is a zero-ton
                continue;
            end
            if floor(Position) == Position
                % The position is an integer, hence it is a singleton
                value = CurrentStream(1,j);
                DFT(Position+1) = value;
                
%% remove the node by subtracting it from all the check nodes it is connected to :
% 
% Consider a stream subsampled at a rate R, the number of elements
% are N/R. The elements in the stream, Xs[0], Xs[1]... can be expressed
% as a sum of the DFT coefficients of the original DFT as below :
%
%   Xs[0] = X[0] + X[N/R] + X[2N/R] + ... + X[N-N/R]
%   Xs[1] = X[1] + X[N/R+1] + X[2N/R+1] + ... + X[N-N/R+1]
%   .
%   .
%   .
%   Xs[N/R-1] = X[N/R-1] + X[N/R+N/R-1] + X[2N/R+N/R-1] + ... + X[N-N/R+N/R-1]
%                                                                   =X[N-1]
% Clearly, every coefficient of the original DFT appears only once in the
% entire set of Xs[0], Xs[1], ... , Xs[N/R-1]. Thus in each stream, a
% particular singleton node solution must be subtracted from 1 particular
% Xs[i] equation in order to remove the solved node

                for k = (1 : len)
                    eqnNo = mod(Position, dftSize/sampleArray(k));
                    % eqnNo tells the location of the positions in a stream
                    % that contain the node at "Position"
                    Stream = FFTsubsampledCell{k};
                    Stream(:,eqnNo+1) = Stream(:,eqnNo+1) - CurrentStream(:,j);
                    FFTsubsampledCell{k} = Stream;
                end
            else
                % It is not a singleton -> multi-ton
            end
        end
    end
    
    finalRun = true;
    for i = (1 : len)
        FFTsubsampledCell{i} = round(FFTsubsampledCell{i},8);
        if any(any(FFTsubsampledCell{i})) ~= 0
            finalRun = false;
        end
    end
    
    if ~finalRun
        DFT = PeelingSolver(FFTsubsampledCell, d, sampleArray, dftSize);
    end
end