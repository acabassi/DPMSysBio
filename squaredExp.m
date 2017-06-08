function covMatrix = squaredExp(logtheta, times)
nRows = length(times);
covMatrix = zeros(nRows, nRows);
se2       = exp(2*logtheta(3));
l2          = exp(2*logtheta(1));
sf2         = exp(2*logtheta(2));
for i = 1:nRows
    for j = 1:nRows
        covMatrix(i,j) = sf2*exp(-(times(i) - times(j))^2/(2*l2));
    end
end
covMatrix = covMatrix  + se2*eye(nRows);
end