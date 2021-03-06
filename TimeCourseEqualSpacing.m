function [output1 output2] = TimeCourseEqualSpacing(input, mode)
switch mode
    case 'init'
        data              = input.data;
        featureNames      = input.featureNames;
        nGenes            = input.nGenes;
        nFeatures         = input.nFeatures;
        timeDiffs         = featureNames - featureNames(1);
        nStartingClusters = ceil(log(nGenes));
        clusterIDs        = random('unid', nStartingClusters, 1, nGenes);
        uniqueIDs         = unique(clusterIDs);
        nStartingClusters = length(uniqueIDs);
        sparseMatrix      = zeros(nGenes,nFeatures);
         % We assume independent normal priors for the log hyperparameters
        
        % Note: order of hypers: length scale, signal variance,
        % noise variance
        hyperPriorParameters = [0, 1; 0, 1; 0, 2]; % [mean s.d.; ...]
        
        % Define the cluster structure
        clusterStruct(1,nStartingClusters+2) = struct(...
            'nGenes', [], ...
            'nFeatures', [], ...
            'timeDiffs', [],...
            'N', [], ...
            'dataCounts', [], ...
            'squaredDataCounts', [], ...
            'logicalGeneIDs', [], ...
            'logMarginalLikelihood',[], ...
            'logHypers', [], ...
            'logPriorOfLogHypers', [], ...
            'squaredHypers', [], ...
            'hyperPriorParams', [], ...
            'firstRowOfCovarianceMatrix', [], ...
            'covarianceMatrixInverses', []);
       
        [clusterStruct.nFeatures       ] = deal(nFeatures);
        [clusterStruct.hyperPriorParams] = deal(hyperPriorParameters);
        [clusterStruct.timeDiffs] = deal(timeDiffs);
        for i = uniqueIDs
            logicalIndices                   = clusterIDs == i;
            indices                          = find(logicalIndices);
            nGenesInCluster                  = length(indices);
            dataInCluster                    = sparseMatrix;
            dataInCluster(indices,:)         = data(logicalIndices,:);
            currentCluster                   = clusterStruct(i);
            currentCluster.logicalGeneIDs    = logicalIndices;  
            currentCluster.dataCounts        = sum(dataInCluster,1);
            currentCluster.squaredDataCounts = sum(dataInCluster.^2,1);
            currentCluster.nGenes            = nGenesInCluster;
            currentCluster.N                 = nFeatures*nGenesInCluster;
            
            logHypers = TimeCourseEqualSpacing(currentCluster, 'sampleHypers');
            currentCluster.logHypers           = logHypers;
            currentCluster.logPriorOfLogHypers = [log(normpdf(logHypers(1)))  log(normpdf(logHypers(2)))  log(normpdf(logHypers(3)))];
            currentCluster.squaredHypers = exp(2*currentCluster.logHypers);
            
            
            l2  = currentCluster.squaredHypers(1);
            sf2 = currentCluster.squaredHypers(2);
            currentCluster.firstRowOfCovarianceMatrix = sf2*exp(-timeDiffs.^2/(2*l2));

            
            currentCluster.covarianceMatrixInverses(1,nGenes) =...
                struct('invertedCovarianceMatrix', [], 'determinant', []);
            
            currentCluster      = TimeCourseEqualSpacing(currentCluster, 'invert');
            currentCluster = TimeCourseEqualSpacing(currentCluster, 'marginal');
            clusterStruct(i) = currentCluster;
        end
        currentCluster = clusterStruct(nStartingClusters+1);
        currentCluster.covarianceMatrixInverses(1,nGenes) =...
            struct('invertedCovarianceMatrix', [], 'determinant', []);
        clusterStruct(nStartingClusters+1) = currentCluster;
        clusterStruct(end) = [];
        
        output1 = clusterStruct;
        output2 = clusterIDs;
    case 'marginal'
        
        nTimes             = input.nFeatures;
        nGenesInCluster    = input.nGenes;
        se2                = input.squaredHypers(3);
        N                  = input.N;
        dataCounter        = input.dataCounts;
        dataSquaredCounter = input.squaredDataCounts;
        L                  = input.covarianceMatrixInverses(nGenesInCluster).invertedCovarianceMatrix;
        invertedMatrix     = reshape(L, nTimes, nTimes);
        logDetK            = input.covarianceMatrixInverses(nGenesInCluster).determinant;
        yColSum            = dataCounter';
        y2                 = sum(dataSquaredCounter);
        input.logMarginalLikelihood = -0.5*(N*log(2*pi)+logDetK+(yColSum'*invertedMatrix*yColSum + y2/se2));
        output1  = input;
        output2  = [];
    case 'invert'
        nGenesInCluster      = input.nGenes;
        se2         = input.squaredHypers(3);%exp(2*input.logHypers(3));
        covFirstRow = input.firstRowOfCovarianceMatrix;
        [invertedK logDetK] =...
            BlockToeplitzInverse(covFirstRow, se2, nGenesInCluster);
        
        if(~isreal(logDetK))  %We should not ever enter here, but just in case
            disp('Numerical error - covariance matrix may not be positive definite')
            logDetK = -realmax;
        end
        input.covarianceMatrixInverses(nGenesInCluster).invertedCovarianceMatrix...
            = invertedK;
        input.covarianceMatrixInverses(nGenesInCluster).determinant = logDetK;
        output1 = input;
        output2  = [];
    case 'sampleHypers'
        logHypers         = zeros(3,1);
        hyperPriorParameters = input.hyperPriorParams;
        for j  = 1:3
            % Initialise the hypers by sampling from the prior
            hyperPriorParams             = hyperPriorParameters(j,:);
            logHypers(j) ...
                = hyperPriorParams(1) + (hyperPriorParams(2)*randn);
        end
        output1 = logHypers;
        output2 = [];
        
    case 'initialiseAuxiliary'
        output1 = input;
        nGenesInCluster    = 1;
        timeDiffs          = output1.timeDiffs;
        output1.nGenes     = nGenesInCluster;
        output1.N          = input.nFeatures;
        
        logHypers             = TimeCourseEqualSpacing(output1, 'sampleHypers');
        output1.logHypers     = logHypers;
        output1.squaredHypers = exp(2*output1.logHypers);
        
        l2                    = output1.squaredHypers(1);
        sf2                   = output1.squaredHypers(2);
        output1.firstRowOfCovarianceMatrix = sf2*exp(-timeDiffs.^2/(2*l2));
        output1      = TimeCourseEqualSpacing(output1, 'invert');
        output2 = [];
end

end