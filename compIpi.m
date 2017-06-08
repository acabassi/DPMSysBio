function sumMat = compIpi(cl, n, PSM)
mat     = cltoSim(cl,n);
sumMat  = sum(sum(tril(mat.*PSM, -1)));