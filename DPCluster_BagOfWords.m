function DPCluster_BagOfWords(fileName, uniqueIdentifier, nSamples, dataType, drawFigures, hyperParameterSamplingFrequency, verbose, initialise, gammaPrior, thinningFreq, inputSeed, GUIflag, varargin)

if(isnan(inputSeed))
    inputSeed = sum(100*clock) + sum(100*double(uniqueIdentifier)); %%clock-seed the random number generator (with a chain-depenedent offset)
end
randn('seed', inputSeed);
rand( 'seed', inputSeed);

timeCourseSwitch   = false;
multinomialSwitch  = false;
bagOfWordsSwitch   = true;

saveFileName = [strtok(fileName, '.'),'_Results_Chain', num2str(uniqueIdentifier)];

if(~initialise)
    % Set up from the provided .mat file:
    savedNSamples                        = nSamples;
    savedDrawFigures                     = drawFigures;
    savedHyperParameterSamplingFrequency = hyperParameterSamplingFrequency;
    savedVerbose                         = verbose;
    load([saveFileName '.mat']);
    nSamples                        = savedNSamples;
    drawFigures                     = savedDrawFigures;
    hyperParameterSamplingFrequency = savedHyperParameterSamplingFrequency;
    verbose                         = savedVerbose;
else
    
    % We need to initialise the variables and structures
    fHandle = @BagOfWords;
    
    % For convenience, store the data in a structure
    allData = importdata(fileName, ',',1);
    data    = allData.data;
    featureNames = allData.textdata(1,2:end);
    geneNames    = allData.textdata(2:end,1);
    
    dataStruct.data         = data;
    dataStruct.geneNames    = geneNames;
    dataStruct.featureNames = featureNames;
    nGenes                  = length(geneNames);
    nFeatures               = length(featureNames);
    dataStruct.nGenes       = length(geneNames);
    dataStruct.nFeatures    = length(featureNames);
    
    % Initialise the clustering partition (structOfClusters)
    [structOfClusters clusterIDs] = feval(fHandle, dataStruct, 'init');
    
    % Set the prior for alpha
    a0 = gammaPrior(1); b0 = gammaPrior(2);  % Parameters for gamma prior
    alpha = gamrnd(a0,1/b0);%1;
    
    
    % We will keep a track of the hyperparameter acceptances
    nHyperProposals = [0 0 0]; nHyperAcceptances = [0 0 0];
    
    % Initialise the output file
    outFile = [saveFileName '.csv'];
    header = 'alpha0,';
    for i = 1:nGenes
        currentGene = geneNames(i);
        header = strcat(header, currentGene, ',');
    end
    header = strcat(header, '\n');
    fid = fopen(outFile, 'wt');
    fprintf(fid, header{1});
    fclose(fid);
    
    % The following is just a variable that enables us to decide how many
    % subplots we need when plotting the clusters
    nClustersOld = 0;
end

nMcmc = nSamples;

if(drawFigures)
    % Draw the initial random configuration
    if(GUIflag == true)
        nClustersOld = doPlots_BagOfWords(clusterIDs, data,featureNames, nClustersOld, 0,varargin{1});
    else
        nClustersOld = doPlots(clusterIDs, data,featureNames, nClustersOld, timeCourseSwitch, multinomialSwitch, bagOfWordsSwitch);
    end
    pause(0.1)
end

for sampleNumber = 1:nMcmc
    if(verbose)
        disp(['Sample number = ', num2str(sampleNumber)]);
    end
    
    % We iterate over the genes:
    for i = 1:dataStruct.nGenes
        % Which clusters are currently occupied?
        occupiedClusterIDs  = unique(clusterIDs);
        nOccupiedClusters   = length(occupiedClusterIDs);
        % Find the cluster in which the current gene resides
        clusterNumber   = clusterIDs(i);
        % Where does this cluster label appear in occupiedClusterIDs ?
        occupiedClusterIndex = find(occupiedClusterIDs == clusterNumber);
        
        % Pick out the data for the current gene
        dataForCurrentGene    = data(i,:);
        
        % Pick out the occupied clusters
        currentClusters   = structOfClusters(occupiedClusterIDs);
        % For each occupied cluster, we will propose adding the gene to
        % each cluster (except for the cluster in which the current gene
        % currently resides, for which we will propose to remove the gene)
        proposedClusters  = currentClusters;
        
        % The following is a switch to indicate if, by removing the current
        % gene from its cluster, we have emptied the cluster
        emptiedClusterSwitch = false;
        for j = (1:nOccupiedClusters)
            % We iterate through the clusters:
            currentClusterLabel = occupiedClusterIDs(j);
            
            currentProposedCluster = proposedClusters(j);
            dataCounts = currentProposedCluster.dataCounts;
            if( currentClusterLabel ~= clusterNumber)
                nGenesInCluster     = currentProposedCluster.nGenes + 1;
                dataCounts          = dataCounts + dataForCurrentGene;
                currentProposedClusterLogicalGeneIDs = currentProposedCluster.logicalGeneIDs;
                currentProposedClusterLogicalGeneIDs(i) = true;
                currentProposedCluster.logicalGeneIDs = currentProposedClusterLogicalGeneIDs;
            else
                nGenesInCluster     = currentProposedCluster.nGenes - 1;
                dataCounts          = dataCounts - dataForCurrentGene;
                currentProposedClusterLogicalGeneIDs = currentProposedCluster.logicalGeneIDs;
                currentProposedClusterLogicalGeneIDs(i) = false;
                currentProposedCluster.logicalGeneIDs = currentProposedClusterLogicalGeneIDs;
            end
            
            currentProposedCluster.nGenes = nGenesInCluster;
            if(nGenesInCluster > 0)
                currentProposedCluster.dataCounts        = dataCounts;
                currentProposedCluster.N                 = nGenesInCluster*nFeatures;
                currentProposedCluster = feval(fHandle, currentProposedCluster, 'marginal'); %TimeCourse(currentProposedCluster, 'marginal'); %getLogMarginalLikelihood(currentProposedCluster);
            else
                emptiedClusterSwitch = true;
            end
            
            proposedClusters(j) = currentProposedCluster;
        end
        
        if(emptiedClusterSwitch)
            % Then we just set the auxiliary component to be the component
            % we just emptied
            auxiliaryCluster = currentClusters(occupiedClusterIndex);
        else
            % We need to initialise the auxiliary component
            auxiliaryCluster      = structOfClusters(end);
            auxiliaryCluster      = feval(fHandle, auxiliaryCluster, 'initialiseAuxiliary');
            auxiliaryCluster.dataCounts        = dataForCurrentGene;
            %%%%%
            auxiliaryClusterLogicalGeneIDs    = false(1, nGenes);
            auxiliaryClusterLogicalGeneIDs(i) = true;
            auxiliaryCluster.logicalGeneIDs = auxiliaryClusterLogicalGeneIDs;
            %%%%%
            
            auxiliaryCluster = feval(fHandle, auxiliaryCluster, 'marginal');
            
            
        end
        logMarginalLikelihoodsWithoutGene = [currentClusters.logMarginalLikelihood];
        logMarginalLikelihoodsWithGene    = [proposedClusters.logMarginalLikelihood];
        
        % We have to swap the entries for the current cluster (since we
        % removed rather than added the gene for this cluster)
        saved = logMarginalLikelihoodsWithoutGene(occupiedClusterIndex);
        logMarginalLikelihoodsWithoutGene(occupiedClusterIndex) = logMarginalLikelihoodsWithGene(occupiedClusterIndex);
        logMarginalLikelihoodsWithGene(occupiedClusterIndex)    = saved;
        
        % Calculate the log marginal likelihood ratios
        logMarginalLikelihoodRatios = logMarginalLikelihoodsWithGene - logMarginalLikelihoodsWithoutGene;
        marginalLikelihoodRatios    = exp(logMarginalLikelihoodRatios);
        
        geneCounts = [currentClusters.nGenes];
        % Again, we have to swap the entries for the current cluster
        geneCounts(occupiedClusterIndex) = proposedClusters(occupiedClusterIndex).nGenes;
        
        % Calculate the conditional posterior probabilities (unnormalised)
        existingClusterProbs  = geneCounts.*marginalLikelihoodRatios;
        newClusterProb        = alpha*exp(auxiliaryCluster.logMarginalLikelihood);
        
        % Normalise these probabilities
        allProbs   = [existingClusterProbs, newClusterProb];
        allProbs   = allProbs/sum(allProbs);
        
        % Determine the cluster to which the gene should be added:
        cumulProbs = cumsum(allProbs);
        
        selected = find(cumulProbs > rand, 1);
        % Update the clusters:
        if(selected~=occupiedClusterIndex)  % We only need to do anything if the gene has been put back in a different cluster
            structOfClusters(clusterNumber) = proposedClusters(occupiedClusterIndex);
            if(selected > nOccupiedClusters)% a new cluster is born...
                % Get a new ID
                available = setdiff((1:(nOccupiedClusters+1)), occupiedClusterIDs);
                updatedClusterLabel = available(1);
                if(updatedClusterLabel == (nOccupiedClusters+1))
                    % Extend the structre array by 1 element
                    structOfClusters(end+1) = structOfClusters(end);
                end
                structOfClusters(updatedClusterLabel) = auxiliaryCluster;
            else
                % Find the appropriate existing ID
                updatedClusterLabel = occupiedClusterIDs(selected);
                structOfClusters(updatedClusterLabel) = proposedClusters(selected);
                
            end
            clusterIDs(i) = updatedClusterLabel;
        end
    end
    
    
    nIndividuals = nGenes;
    %%% Having completed the Gibbs sampling for the component indicator
    %%% variables, we now sample a new alpha (see Escobar and West, 1995)
    k    = length(unique(clusterIDs));
    eta  = betarnd((alpha+1), nIndividuals);
    
    %Escobar & West gives (pi_eta/(1 - pi_eta)) = A/B, where A & B are:
    A    = a0 + k - 1;
    B    = nIndividuals*(b0 - log(eta));
    
    %Rearranging, we have pi_eta = A/(A+B)
    pi_eta     = A/(A + B);
    
    if rand < pi_eta
        alpha = gamrnd(a0+k  , 1/(b0 - log(eta)));
    else
        alpha = gamrnd(a0+k-1, 1/(b0 - log(eta)));
    end
    
    
    %%% save down the current sample (only save every 'thinningFreq'-th
    %%% sample)
    if(mod(sampleNumber,thinningFreq) == 0)
        occupiedClusterIDs  = unique(clusterIDs);
        outputVector = [alpha, clusterIDs];
        dlmwrite(outFile, outputVector, '-append', 'delimiter', ',');
    end
    
    
    if(mod(sampleNumber,thinningFreq) == 0 || sampleNumber == nMcmc)
        clear savedDrawFigures savedHyperParameterSamplingFrequency savedNSamples savedVerbose;
        save([saveFileName '.mat']);
    end
    if(verbose)
        disp(['alpha = ', num2str(alpha)]);
        disp(['nClusters =', num2str(length(unique(clusterIDs)))])
    end
    
    if(drawFigures)
        refresh
        if(GUIflag == true)
            nClustersOld = doPlots_BagOfWords(clusterIDs, data,featureNames, nClustersOld, sampleNumber, varargin{1});
        else
            nClustersOld = doPlots(clusterIDs, data,featureNames, nClustersOld, timeCourseSwitch, multinomialSwitch, bagOfWordsSwitch);
        end
        pause(0.1)
    end
    
    
end

end