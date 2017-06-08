function [invertedBlockMatrix blockMatrixDeterminant] = invertBlockMatrix(covarianceMatrix, nGenes, nTimes, noise)

% Assume that the elements of the covariance matrix have *not* been ordered
% by time (so that indices 1:nTimes correspond to Gene 1, then the next
% nTimes rows/columns correspond to Gene 2, ... and so on).
%noiselessDiagonal = covarianceMatrix((nTimes+1),1);
covarianceMatrix2 = covarianceMatrix((nTimes+1):(2*nTimes), 1:nTimes);
% We first want to "collapse" the matrix, to avoid unnecessary storage
lowerDiagonal = covarianceMatrix2(logical(tril(ones(nTimes))));
% The below allows us to reconstruct the covariance matrix
collapsedCovarianceMatrix = [lowerDiagonal; noise];


[invertedCollapsedBlockMatrix blockMatrixDeterminant] = invertCollapsedBlockMatrix(collapsedCovarianceMatrix, nGenes, nTimes,noise);

% Reconstruct the ``uncollapsed" matrix

matrix1 = zeros(nTimes, nTimes);
i = logical(tril(ones(nTimes, nTimes)));
matrix1(i) = invertedCollapsedBlockMatrix;
matrix1 = matrix1.';
matrix1(i) = invertedCollapsedBlockMatrix;


invertedBlockMatrix =  repmat(matrix1,nGenes,nGenes) + eye(nGenes*nTimes)/noise;
%
% USE THE IDENTITY FOR SUB-DIVIDING A BLOCK MATRIX
% K = [A B]    and E = D - C*A_inverse*B
%     [C D]


%nTimes = size(covarianceMatrix,1); % Since the matrix is square, this is also nColumns
%lowerDiagonal = covarianceMatrix(tril(covarianceMatrix,-1)>0);
%diagonalElements  = covarianceMatrix(1,1);
%invertedBlockMatrix = invertBlockMatrixVectorRepresentation(lowerDiagonal, diagonalElements, nGenes, nRows);
%A = covarianceMatrix(1:nGenes, 1:nGenes)

end