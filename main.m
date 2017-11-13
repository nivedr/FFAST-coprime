
% Input Sparse DFT array with:
% k = 3
% N1 = 336 = 16*21
% N2 = 323 = 17*19
% Rate of subsampling: 16 (vertically) * 17 (horizontally) -> /272
%                      21 (vertically) * 19 (horizontally) -> /399

inputArray = zeros(304,306);
inputArray(3,3) = 10;
inputArray(233,100) = 23;
inputArray(145,301) = 1;
inputArray(33,166) = 1.5;

% Performing DFT on the iDFT of the sparse array. (Implicit iDFT performed)
DFT = FFAST(inputArray,304,306,2,[16,19],[18,17]);

working = any(any((inputArray - DFT)));