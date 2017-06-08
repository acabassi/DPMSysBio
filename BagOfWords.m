function [output1 output2] = BagOfWords(input, mode)
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
            'hyperParameters', []);
       
        [clusterStruct.nFeatures       ] = deal(nFeatures);
        hyperParameters = 0.5 + zeros(1,nFeatures);
        
        [clusterStruct.hyperParameters ] = deal(hyperParameters);
        
        
        for i = uniqueIDs
            logicalIndices                   = clusterIDs == i;
            indices                          = find(logicalIndices);
            nGenesInCluster                  = length(indices);
            dataInCluster                    = sparseMatrix;
            dataInCluster(indices,:)         = data(logicalIndices,:);
            currentCluster                   = clusterStruct(i);
            currentCluster.logicalGeneIDs    = logicalIndices;  
            currentCluster.dataCounts        = sum(dataInCluster);
            currentCluster.nGenes            = nGenesInCluster;
            currentCluster.N                 = nFeatures*nGenesInCluster;
            currentCluster                   = BagOfWords(currentCluster, 'marginal');
            clusterStruct(i) = currentCluster;
        end
        currentCluster = clusterStruct(nStartingClusters+1);
        clusterStruct(nStartingClusters+1) = currentCluster;
        clusterStruct(end) = [];
        
        output1 = clusterStruct;
        output2 = clusterIDs;
    case 'marginal'
        dataCounts   = input.dataCounts;
        nTotalCounts = sum(dataCounts);
        beta         = input.hyperParameters;
        sumBeta      = sum(beta);
        
        logMarginalLikelihood = sum(gammaln(dataCounts+beta));
        logMarginalLikelihood = logMarginalLikelihood + gammaln(sumBeta);
        logMarginalLikelihood = logMarginalLikelihood - sum(gammaln(beta));
        logMarginalLikelihood = logMarginalLikelihood - gammaln(sumBeta + nTotalCounts);

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







