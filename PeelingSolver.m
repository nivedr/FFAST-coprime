%% Function PeelingSolver
% Solves for the DFT coefficients of the original DFT from the subsampled
% DFT coefficients using the Peeling Decoder. This function is used in the
% <FFAST.html 2d-FFAST> implmentation.

%% Parameters of the function
%
% <html>
% <p><pre><h3>Inputs</h3>
% Ystream       : coefficients of subsampled DFT
% N1            : Input column dimension size
% N2            : Input row dimension size
% N1sampleArray : Subsampling rates along Y
% N1sampleArray : Subsampling rates along X<br>
% <h3>Outputs</h3>
% DFT           : Computed DFT (from subsampled versions)
% </pre></p></html> 

function [DFT] = PeelingSolver(YStream, N1, N2, N1sampleArray, N2sampleArray)

    
%% Implementation
% The implementation of the Peeling Decoder has been illustrated below.
% This function is used in conjuction with the  <FFAST.html 2d-FFAST>
%
%%
% Check to see if the FFT of the subsampled array has been set/reduced to
% all zeros. This is a check that gives whether the iterations have succesfully
% solved for the DFT coefficients of the original 2d signal. 

    % Initialize output
    DFT = zeros(N1, N2);

    % Check if subsampled streams are empty using the "empty" flag
    empty = true; 
    for i = (1 : length(N1sampleArray))
        cell = YStream{i};
        for p = (1 : N1/N1sampleArray(i))
            for q = (1 : N2/N2sampleArray(i))
                cellContent = cell{p,q};
                if any(round(cellContent, 8)) ~= 0
                    empty = false;
                end
            end
        end
    end

    if empty
        return;
    end

%% Peeling Decoder

    % Iterating over streams
    for i = (1 : length(N1sampleArray))
        P = YStream{i};     % 21*19 cell each of which holds a 3*1 vector
        
        % Iterating over columns in a stream
        for j = (1 : N1/N1sampleArray(i))
            
            % Iterating over rows in a stream
            for k = (1 : N2/N2sampleArray(i))
                
                % Pick up every "support" in the stream that corresponds to 3 elements:
                % Sum of a set of DFT coefficients (of the original DFT)
                % Sum of phase shifted versions of the DFT coefficients (along Y)
                % Sum of phase shifted versions of the DFT coefficients (along X)
                support = round(P{j,k},8); % round for numerical stability
                
                % Check if the node has value 0 (In which case it CANNOT be
                % a singleton)
                if support(1,1) ~= 0
                    % Singleton Test :
                    % compute estimated row & column. If both are whole
                    % numbers, high probability that it is a singleton node
                    % round for numerical stability
                    column = round((N1/(2*pi))*angle(support(2,1)*conj(support(1,1))), 4);
                    row = round((N2/(2*pi))*angle(support(3,1)*conj(support(1,1))), 4);
                else
                    % The node is a zeroton or a multiton
                    continue;
                end
                
                % Compute the value of the node and remove it and its phase
                % shifted versions from the supports in which it occurs
                % ** if it is a singleton
                if floor(column) == column && floor(row) == row % singleton test
                    value = support(1,1);
                    % Update output
                    DFT(mod(column,N1)+1, mod(row,N2)+1) = value;
                
%% remove the node by subtracting it from all the check nodes it is connected to :
% 
% Consider a stream subsampled at a rate (Ry,Rx), the number of elements
% are $N1  Ry, N2  Rx$. The elements in the stream, Xs[00,10,01], Xs[00,10,01]
% (00->regular, 10->Y-delayed(Y-phase shifted), 01->X-delayed(X-phase shifted))
% can be expressed as a sum of the DFT coefficients of the original DFT as 
% below :
%
%   Xs00[yx] =     X[y,x]   +    X[N1/Ry+y,x]    +    X[2N1/Ry+y,x]    + ... + X[N1-N1/Ry+y,x]
%            + X[y,N2/Rx+x] + X[N1/Ry+y,N2/Rx+x] + X[2N1/Ry+y,N2/Rx+x] + ... + X[N1-N1/Ry+y,N2/Rx+x]
%            ...
% (y , x) \in (0 : N1/Ry-1 , 0 : N2/Rx-1)
% Clearly, every coefficient of the original DFT appears only once in the
% entire set of Xs[0], Xs[1], ... , Xs[N/R-1] of a stream
% The y,x where it occurs can be found as (k1.N1/Ry+y, k2.N2/Rx+x) = (col,row)
% or : y = col|N1/Ry, x = col|N2/Rx
%
                    for l = (1 : length(N1sampleArray))
                        colNo = mod(column, N1/N1sampleArray(l));
                        rowNo = mod(row, N2/N2sampleArray(l));
                        % remove the singleton solution
                        P1 = YStream{l};
                        P1{colNo+1, rowNo+1} = P1{colNo+1, rowNo+1} - support;
                        % update the subsampled FFT array back
                        YStream{l} = P1;
                    end
                end    
            end
        end
    end

    finalRun = true;
    for i = (1 : length(N1sampleArray))
        cell = YStream{i};
        for p = (1 : N1/N1sampleArray(i))
            for q = (1 : N2/N2sampleArray(i))
                cellContent = cell{p,q};
                if any(round(cellContent, 8)) ~= 0
                    finalRun = false;
                end
            end
        end
    end
        
    if ~finalRun
        DFT = PeelingSolver(YStream, N1, N2, N1sampleArray, N2sampleArray);
    end
 end   
