function [output1 output2] = Multinomial(input, mode)
switch mode
    case 'init'
        data              = input.data;
        nGenes            = input.nGenes;
        nFeatures         = input.nFeatures;
        
        nStartingClusters = ceil(log(nGenes));
        clusterIDs        = random('unid', nStartingClusters, 1, nGenes);
        uniqueIDs         = unique(clusterIDs);
        nStartingClusters = length(uniqueIDs);
        sparseMatrix      = zeros(nGenes,nFeatures);
 
        % Define the cluster structure
        clusterStruct(1,nStartingClusters+2) = struct(...
            'nGenes', [], ...
            'nFeatures', [], ...
            'N', [], ...
            'dataCounts', [], ...
            'logicalGeneIDs', [], ...
            'logMarginalLikelihood', [],...
            'hyperParameters', [],...
            'dataLevels', [],...
            'dataCountIndexHelper', []);
       
        [clusterStruct.nFeatures       ] = deal(nFeatures);
        dataLevels      = unique(data);
        dataLevels      = reshape(dataLevels, 1, length(dataLevels)); %ensure we have a row vector
        dataProportions = histc(data(:), dataLevels);
        dataProportions = dataProportions/sum(dataProportions);
        hyperParameters = length(dataProportions)*0.5*dataProportions;
        
        [clusterStruct.hyperParameters ] = deal(hyperParameters);
        nLevels = length(dataLevels);
        dataCountIndexHelper = 0:nLevels:(nLevels*(nFeatures-1));
        [clusterStruct.dataCountIndexHelper ] = deal(dataCountIndexHelper);
        [clusterStruct.dataLevels] = deal(dataLevels);
        for i = uniqueIDs
            logicalIndices                   = clusterIDs == i;
            indices                          = find(logicalIndices);
            nGenesInCluster                  = length(indices);
            dataInCluster                    = sparseMatrix;
            dataInCluster(indices,:)         = data(logicalIndices,:);
            currentCluster                   = clusterStruct(i);
            %currentCluster.data              = dataInCluster;
            currentCluster.logicalGeneIDs    = logicalIndices;  
            currentCluster.dataCounts        = histc(dataInCluster(indices,:),dataLevels);
            currentCluster.nGenes            = nGenesInCluster;
            currentCluster.N                 = nFeatures*nGenesInCluster;
            currentCluster                   = Multinomial(currentCluster, 'marginal');
            clusterStruct(i) = currentCluster;
        end
        currentCluster = clusterStruct(nStartingClusters+1);
        clusterStruct(nStartingClusters+1) = currentCluster;
        clusterStruct(end) = [];
        
        output1 = clusterStruct;
        output2 = clusterIDs;
    case 'marginal'
        hyperParameters = input.hyperParameters;
        nFeatures       = input.nFeatures;
        dataCounts     = input.dataCounts;
        nGenes         = input.nGenes;
        beta      = repmat(hyperParameters, 1, nFeatures);
        sumBeta   = sum(beta, 1);
        logMarginalLikelihood = sum(sum(gammaln(dataCounts+beta)));
        logMarginalLikelihood = logMarginalLikelihood + sum(gammaln(sumBeta));
        logMarginalLikelihood = logMarginalLikelihood - sum(sum(gammaln(beta)));
        logMarginalLikelihood = logMarginalLikelihood - sum(gammaln(nGenes + sumBeta));
        input.logMarginalLikelihood = logMarginalLikelihood;
        output1 = input;
    case 'initialiseAuxiliary'
        output1 = input;
        nGenesInCluster    = 1;
        output1.nGenes     = nGenesInCluster;
        output1.N          = input.nFeatures;
        output2 = [];
end

end







