function [y logDet] = BlockToeplitzInverse(r, noise, blockSize)
% Here, r must be the first row of the covariance matrix that we 
% wish to invert (in "collapsed" form).

n = length(r);
BlockOnDiagonal = r(1);
[InverseBlockOnDiagonal newNoise logDetBlockOnDiagonal] = blockInverse(BlockOnDiagonal, noise, blockSize);


rOld = r;
r(1) = [];

r = (blockSize*InverseBlockOnDiagonal + newNoise)*r;
% Call BlockTrench, which assumes that the matrix on the diagonal is the
% identity
[y logDet] = BlockTrench(r, blockSize);


y = (blockSize*InverseBlockOnDiagonal + newNoise)*y;
correctionTerm = sparse(zeros(1,n));
correctionTerm(1:n) = InverseBlockOnDiagonal;
y = y + diag(correctionTerm);
logDet = (n*logDetBlockOnDiagonal) + logDet;
