function myMat = cltoSim(cl, n)
[X, Y] = meshgrid(cl,cl);
X = X(:); Y = Y(:);
myMat = reshape(X==Y, n, n);
end

