addpath('/Users/paulk/MatlabCode/DPCluster/gpml-matlab/gpml-matlab/gpml/');
N = 17*output(1).nGenes;
cluster = output(1);
inputs = repmat([0:10:160]', cluster.nGenes, 1);

covarianceMatrix = zeros(N,N);

noise       = exp(2*cluster.logHypers(3));
timeDiffs   = cluster.timeDiffs;
l2          = exp(2*cluster.logHypers(1));
sf2         = exp(2*cluster.logHypers(2));

for i = 1:N
    for j = 1:N
        covarianceMatrix(i,j) = inputs(i) - inputs(j);
    end
end

covarianceMatrix = covarianceMatrix.^2;
covarianceMatrix = -covarianceMatrix./(2*l2);
covarianceMatrix = sf2.*exp(covarianceMatrix);
covarianceMatrix = covarianceMatrix + noise*eye(N);
inverted = inv(covarianceMatrix);


tdata    = cluster.data';
outputs  = tdata(:);
logtheta = cluster.logHypers;
covfunc  = {'covSum', {'covSEiso', 'covNoise'}};
x = inputs; y = outputs;
out1 = gpr(logtheta, covfunc, x, y);





nGenes      = cluster.nGenes;
noise       = exp(2*cluster.logHypers(3));
timeDiffs   = cluster.timeDiffs;
l2          = exp(2*cluster.logHypers(1));
sf2         = exp(2*cluster.logHypers(2));
covFirstRow = sf2*exp(-timeDiffs.^2/(2*l2));
[invertedK logDetK] =...
    BlockToeplitzInverse(covFirstRow, noise, nGenes);

%Note:
%inv(covarianceMatrix) = invertedK + (eye(N)./noise)


covM = zeros(17,17);
for i = 1:17
    for j = 1:17
        td = inputs(i)-inputs(j);
        covM(i,j) = sf2*exp(-(td^2)/(2*l2));
    end
end
covM = covM + noise*eye(17);

[invertedK logDetK] =...
    BlockToeplitzInverse(covFirstRow, noise, 1);

logDetK
%0.667479210412824
log(det(covM))
%0.667479210412823

%inv(covM) == invertedK + (eye(17)./noise);



