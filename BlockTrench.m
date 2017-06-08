function [B logDet ] = BlockTrench(r,blockSize)
n = 1 + length(r);
B = sparse(zeros(n,n));
[y logDet] = BlockDurbin(r,blockSize);
reciprocalBlockSize = 1/blockSize;
ry = r*y(1:(n-1))*blockSize;
gamma =  (-ry/(ry*blockSize + 1));
gammaPlusReciprocalBlockSize = gamma + reciprocalBlockSize;

v = ((gamma*blockSize+1)*y((n-1):-1:1));
B(1,1) = gamma;
B(1,2:n) = v((n-1):-1:1)';
%B(2:n,1) = B(1,2:n);
%B(n,2:n) = B(1,(n-1):-1:1);
%B(2:(n-1),n) = B(n, 2:(n-1)); 
floored = floor((n-1)/2);
for i = 2:(1+floored)
    for j = i:(n-i+1)
        B(i,j) = B(i-1,j-1) + (   (v(n+1-j)*v(n+1-i)) - (v(i-1)*v(j-1)))/gammaPlusReciprocalBlockSize;
%        B(j,i) = B(i,j);
%        B(n-j+1, n-i+1) = B(i,j);
%        B(n-i+1, n-j+1) = B(n-j+1, n-i+1);
    end
end

BT = B';
B = B + BT;
B = B + rot90(B,2);
B = full(B);
B(1:n+1:end) = B(1:n+1:end)/2;
B(n:n-1:end-1) = B(n:n-1:end-1)/2;