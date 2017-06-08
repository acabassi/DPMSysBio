function [y l] = BlockDurbin(r, blockSize)
% Note: we assume identity matrices down the diagonal
n = length(r);
y = zeros(n,1);
y(1) = -r(1); beta = 1; alpha = -r(1);
l = 0;
for k = 1:(n-1)
    beta = beta - beta*alpha*alpha*blockSize*blockSize;
    l = l + log(beta);
    %alpha = -(r(k+1) + blockDotProduct(r(k:-1:1), y(1:k), blockSize))/beta;
    alpha = -(r(k+1) +  r(k:-1:1)*y(1:k)*blockSize   )/beta;
    z = y(1:k) + blockSize*alpha*y(k:-1:1);
    y(1:k+1) = [z;alpha];
end
beta = beta - beta*alpha*alpha*blockSize*blockSize;
l = l + log(beta);