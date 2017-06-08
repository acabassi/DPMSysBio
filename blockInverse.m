function [y newnoise logDet] = blockInverse(r, noise, blockSize)
newnoise = 1/noise;
y = newnoise*(-r/(r*blockSize + noise));
% detr = (noise^(blockSize-1))*((blockSize*r)+noise)
logDet = (blockSize-1)*log(noise) + log((blockSize*r)+noise);
% Therefore, if blockSize*r = -noise, we are in trouble!  [This seems
% unlikely to happen in practice]