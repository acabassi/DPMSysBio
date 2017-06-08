function DPCluster(fileName, uniqueIdentifier, nSamples, dataType, drawFigures, hyperParameterSamplingFrequency, verbose, initialise, inputSeed)

if(isnan(inputSeed))
    inputSeed = sum(100*clock) + sum(100*double(uniqueIdentifier));
    randn('seed', inputSeed);%%clock-seed the random number generator (with a chain-depenedent offset)
    rand('seed', inputSeed);%%clock-seed the random number generator (with a chain-depenedent offset)
end

saveFileName = [strtok(fileName, '.'),'_Results_Chain', num2str(uniqueIdentifier)];

if(~initialise)
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
        
    % 1. Specify the data-type specific function to initialise the struct
    switch dataType
        case 'TimeCourseEqualSpacing'
            fHandle = @TimeCourse;
        case 'TimeCourseUnequalSpacing'
            fHandle = @TimeCourse;
        case 'TimeCourse'
            fHandle = @TimeCourse;
        case 'DiscreteGeneExpression'
            fHandle = @Multinomial;
        case 'Multinomial'
            fHandle = @Multinomial;
        case 'ChIPchip'
            fHandle = @BagOfWords;
        case 'TranscriptionFactor'
            fHandle = @BagOfWords;
        case 'BagOfWords'
            fHandle = @BagOfWords;
    end
    
    timeCourseSwitch   = logical(sum(strcmp(dataType, {'TimeCourseEqualSpacing','TimeCourse' , 'TimeCourseUnequalSpacing'})));
    multinomialSwitch  = logical(sum(strcmp(dataType, {'Multinomial','DiscreteGeneExpression'})));
    bagOfWordsSwitch   = logical(sum(strcmp(dataType, {'BagOfWords','TranscriptionFactor','ChIPchip'})));
    
    
    % 2. Read in data (genenames, featurenames, and the data)
    allData = importdata(fileName, ',',1);
    data    = allData.data;
    featureNames = allData.textdata(1,2:end);
    geneNames    = allData.textdata(2:end,1);
    
    if(timeCourseSwitch)
        featureNames = cellfun(@str2num,featureNames);
    end
    
    if(multinomialSwitch)
        %We require the data to be numbers 1,2, ..., nLevels
        dataLevels = unique(data)';
        nLevels = length(dataLevels);
        if(~isequal(dataLevels, 1:nLevels))
            %Force the required format
            data = data - max(data(:));
            dataLevels = unique(data)';
            dataLevels = dataLevels(end:-1:1);
            requiredLevels = nLevels:-1:1;
            for i = 1:nLevels
                currentLevel = dataLevels(i);
                requiredLevel = requiredLevels(i);
                data(data == currentLevel) = requiredLevel;
            end
            
        end
    end
        
    dataStruct.data         = data;
    dataStruct.geneNames    = geneNames;
    dataStruct.featureNames = featureNames;
    nGenes                  = length(geneNames);
    nFeatures               = length(featureNames);
    dataStruct.nGenes       = length(geneNames);
    dataStruct.nFeatures    = length(featureNames);
    
    
    [structOfClusters clusterIDs] = feval(fHandle, dataStruct, 'init');
    
    a0 = 2; b0 = 4;  % Parameters for gamma prior
    alpha = gamrnd(a0,1/b0);%1;

    sparseMatrix      = zeros(nGenes,nFeatures);
    
    
    % 3. Gibbs sample
    nHyperProposals = [0 0 0]; nHyperAcceptances = [0 0 0];
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
    
    nClustersOld = 0;
end

nMcmc = nSamples;

if(drawFigures)
    figure
    nClustersOld = doPlots(clusterIDs, data,featureNames, nClustersOld, timeCourseSwitch, multinomialSwitch, bagOfWordsSwitch);
    pause(0.1)
end

for sampleNumber = 1:nMcmc
    if(verbose)
%         tic
        disp(['Sample number = ', num2str(sampleNumber)]);
    end

    for i = 1:dataStruct.nGenes
        % Which clusters are currently occupied?
        occupiedClusterIDs  = unique(clusterIDs);
        nOccupiedClusters   = length(occupiedClusterIDs);
        % Find the cluster in which the current gene resides
        clusterNumber   = clusterIDs(i);
        % Where does this cluster label appear in occupiedClusterIDs ?
        occupiedClusterIndex = find(occupiedClusterIDs == clusterNumber);
        
        % Pick out the appropriate cluster structure
        currentCluster    = structOfClusters(clusterNumber);
        dataForCurrentGene    = data(i,:);
        
        currentClusters   = structOfClusters(occupiedClusterIDs);
        proposedClusters  = currentClusters;
        
        emptiedClusterSwitch = false;
        for j = (1:nOccupiedClusters)
            currentClusterLabel = occupiedClusterIDs(j);
            
            currentProposedCluster = proposedClusters(j);
            %%%%%%%%%%%%%%%%%%%%%
            if(timeCourseSwitch)
                if( currentClusterLabel ~= clusterNumber)
                    nGenesInCluster   = currentProposedCluster.nGenes + 1;
                    dataCounts        = currentProposedCluster.dataCounts + dataForCurrentGene;
                    squaredDataCounts = currentProposedCluster.squaredDataCounts+ dataForCurrentGene.^2;
                    currentProposedClusterLogicalGeneIDs = currentProposedCluster.logicalGeneIDs;
                    currentProposedClusterLogicalGeneIDs(i) = true;
                    currentProposedCluster.logicalGeneIDs = currentProposedClusterLogicalGeneIDs;
                else
                    nGenesInCluster   = currentProposedCluster.nGenes - 1;
                    dataCounts        = currentProposedCluster.dataCounts - dataForCurrentGene;
                    squaredDataCounts = currentProposedCluster.squaredDataCounts- dataForCurrentGene.^2;
                    currentProposedClusterLogicalGeneIDs = currentProposedCluster.logicalGeneIDs;
                    currentProposedClusterLogicalGeneIDs(i) = false;
                    currentProposedCluster.logicalGeneIDs = currentProposedClusterLogicalGeneIDs;
                end
                
                currentProposedCluster.nGenes            = nGenesInCluster;
                if(nGenesInCluster > 0)
                    currentProposedCluster.dataCounts        = dataCounts;
                    currentProposedCluster.squaredDataCounts = squaredDataCounts;
                    currentProposedCluster.N                 = nGenesInCluster*nFeatures;
                    if(isempty(currentProposedCluster.covarianceMatrixInverses(nGenesInCluster).determinant))
                        %[invertedK logDetK] = invertCovarianceMatrix(currentProposedCluster);
                        currentProposedCluster = feval(fHandle, currentProposedCluster, 'invert');
                        % We also want to update the covariance matrices of the
                        % currentClusters too, in order to avoid recalculation
                        structOfClusters(currentClusterLabel).covarianceMatrixInverses(nGenesInCluster) = currentProposedCluster.covarianceMatrixInverses(nGenesInCluster);
                    end
                    currentProposedCluster = feval(fHandle, currentProposedCluster, 'marginal'); %TimeCourse(currentProposedCluster, 'marginal'); %getLogMarginalLikelihood(currentProposedCluster);
                else
                    emptiedClusterSwitch = true;
                end
                
            end
            %%%%%%%%%%%%%%%%%%%%%
            if(multinomialSwitch)
                dataCountIndexHelper = currentProposedCluster.dataCountIndexHelper;
                dataCounts = currentProposedCluster.dataCounts;
                indices = dataCountIndexHelper+dataForCurrentGene;
                if( currentClusterLabel ~= clusterNumber)
                    nGenesInCluster   = currentProposedCluster.nGenes + 1;
                    dataCounts(indices) = dataCounts(indices) + 1;
                    currentProposedClusterLogicalGeneIDs = currentProposedCluster.logicalGeneIDs;
                    currentProposedClusterLogicalGeneIDs(i) = true;
                    currentProposedCluster.logicalGeneIDs = currentProposedClusterLogicalGeneIDs;
                else
                    nGenesInCluster   = currentProposedCluster.nGenes - 1;
                    dataCounts(indices) = dataCounts(indices) - 1;
                    currentProposedClusterLogicalGeneIDs = currentProposedCluster.logicalGeneIDs;
                    currentProposedClusterLogicalGeneIDs(i) = false;
                    currentProposedCluster.logicalGeneIDs = currentProposedClusterLogicalGeneIDs;
                end
                
                currentProposedCluster.nGenes            = nGenesInCluster;
                if(nGenesInCluster > 0)
                    currentProposedCluster.dataCounts        = dataCounts;
                    currentProposedCluster.N                 = nGenesInCluster*nFeatures;
                    currentProposedCluster = feval(fHandle, currentProposedCluster, 'marginal'); %TimeCourse(currentProposedCluster, 'marginal'); %getLogMarginalLikelihood(currentProposedCluster);
                else
                    emptiedClusterSwitch = true;
                end
            end
            %%%%%%%%%%%%%%%%%%%%%
            if(bagOfWordsSwitch)
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
            end
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            
            
            proposedClusters(j) = currentProposedCluster;
        end
        
        if(emptiedClusterSwitch)
            auxiliaryCluster = currentClusters(occupiedClusterIndex);
        else
            auxiliaryCluster      = structOfClusters(end);
            auxiliaryCluster      = feval(fHandle, auxiliaryCluster, 'initialiseAuxiliary');
            if(timeCourseSwitch)
                auxiliaryCluster.dataCounts        = dataForCurrentGene;
                auxiliaryCluster.squaredDataCounts = dataForCurrentGene.^2;
            end
            if(multinomialSwitch)
                auxiliaryCluster.dataCounts        = histc(dataForCurrentGene, auxiliaryCluster.dataLevels,1);
            end
            if(bagOfWordsSwitch)
                auxiliaryCluster.dataCounts        = dataForCurrentGene;
            end
            %%%%%
            auxiliaryClusterLogicalGeneIDs    = false(1, nGenes);
            auxiliaryClusterLogicalGeneIDs(i) = true;
            auxiliaryCluster.logicalGeneIDs = auxiliaryClusterLogicalGeneIDs;
            %%%%%
            
            auxiliaryCluster = feval(fHandle, auxiliaryCluster, 'marginal');
            
            
        end
        logMarginalLikelihoodsWithoutGene = [currentClusters.logMarginalLikelihood];
        logMarginalLikelihoodsWithGene    = [proposedClusters.logMarginalLikelihood];
        
        % We have to swap the entries for the current cluster:
        saved = logMarginalLikelihoodsWithoutGene(occupiedClusterIndex);
        logMarginalLikelihoodsWithoutGene(occupiedClusterIndex) = logMarginalLikelihoodsWithGene(occupiedClusterIndex);
        logMarginalLikelihoodsWithGene(occupiedClusterIndex)    = saved;
        
        logMarginalLikelihoodRatios = logMarginalLikelihoodsWithGene - logMarginalLikelihoodsWithoutGene;
        marginalLikelihoodRatios    = exp(logMarginalLikelihoodRatios);
        
        geneCounts = [currentClusters.nGenes];
        % Again, we have to swap the entries for the current cluster
        geneCounts(occupiedClusterIndex) = proposedClusters(occupiedClusterIndex).nGenes;
        
        existingClusterProbs  = geneCounts.*marginalLikelihoodRatios;
        newClusterProb        = alpha*exp(auxiliaryCluster.logMarginalLikelihood);
        
        allProbs   = [existingClusterProbs, newClusterProb];
        allProbs   = allProbs/sum(allProbs);
        
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
    
    
    
    %%% sample alpha
    k = length(unique(clusterIDs));
    eta = betarnd((alpha+1), nGenes);
    
    if rand < eta
        alpha = gamrnd(a0+k, 1/(b0 - log(eta)));
    else
        alpha = gamrnd(a0+k-1, 1/(b0 - log(eta)));
    end
    
    
    %%% save down the current sample (only save every 5th sample)
    if(mod(sampleNumber,5) == 0)
        occupiedClusterIDs  = unique(clusterIDs);
        if(timeCourseSwitch)
            currentLogHypers = [structOfClusters(occupiedClusterIDs).logHypers];
            currentLogHypers = currentLogHypers(:);
            outputVector = [alpha, clusterIDs, currentLogHypers'];
        else
            outputVector = [alpha, clusterIDs];
        end
        dlmwrite(outFile, outputVector, '-append', 'delimiter', ',');
    end
    
    %%% sample hypers
    if(timeCourseSwitch)
        if(mod(sampleNumber,hyperParameterSamplingFrequency) == 0)
            % Which clusters are currently occupied?
            occupiedClusterIDs  = unique(clusterIDs);
            nOccupiedClusters   = length(occupiedClusterIDs);
            
            emptyCovMatStruct = structOfClusters(end).covarianceMatrixInverses;
            
            for i = 1:nOccupiedClusters
                
                %%MODIFICATION: update each of the hypers in turn (rather than
                %%all at once)
                for j = 1:3
                    nHyperProposals(j) = nHyperProposals(j) + 1;
                    clusterNumber    = occupiedClusterIDs(i);
                    currentCluster   = structOfClusters(clusterNumber);
                    currentLogHypers = currentCluster.logHypers;
                    currentlogPriorOfLogHypers = currentCluster.logPriorOfLogHypers;
                    currentLogMarginalLikelihood = currentCluster.logMarginalLikelihood;
                    
                    proposedCluster   = currentCluster;
                    proposedCluster.covarianceMatrixInverses = emptyCovMatStruct;
                    proposedLogHypers = currentLogHypers;
                    if (j == 1)
                        proposedLogHypers(j) = proposedLogHypers(j) + (1.0*randn);
                    else
                        proposedLogHypers(j) = proposedLogHypers(j) + (0.5*randn);
                    end
                    proposedCluster.logHypers = proposedLogHypers;
                    squaredHypers = currentCluster.squaredHypers;
                    squaredHypers(j) = exp(2*proposedLogHypers(j));
                    proposedCluster.squaredHypers = squaredHypers;
                    l2          = proposedCluster.squaredHypers(1);
                    sf2         = proposedCluster.squaredHypers(2);
                    timeDiffs   = proposedCluster.timeDiffs;
                    if(strcmp(dataType, 'TimeCourseEqualSpacing'))
                        proposedCluster.firstRowOfCovarianceMatrix = sf2*exp(-timeDiffs.^2/(2*l2));
                    else
                        lowerTriangularLogicalMatrix = proposedCluster.lowerTriangularLogicalMatrix;
                        covarianceMatrix = sf2*exp(timeDiffs/(2*l2));
                        lowerTriangularPart = covarianceMatrix(lowerTriangularLogicalMatrix);
                        proposedCluster.lowerTriangularPartOfCovarianceMatrix = lowerTriangularPart;
                    end
                    
                    proposedLogPriorOfLogHypers = currentlogPriorOfLogHypers;
                    proposedLogPriorOfLogHypers(j) = log(normpdf(proposedLogHypers(j)));
                    proposedCluster.logPriorOfLogHypers = proposedLogPriorOfLogHypers;
                    proposedCluster = feval(fHandle, proposedCluster, 'invert');
                    proposedCluster = feval(fHandle, proposedCluster, 'marginal');
                    proposedLogMarginalLikelihood = proposedCluster.logMarginalLikelihood;
                    
                    logRatio = proposedLogMarginalLikelihood + ...
                        sum(proposedLogPriorOfLogHypers) - ...
                        currentLogMarginalLikelihood - ...
                        sum(currentlogPriorOfLogHypers);
                    
                    if rand < exp(logRatio)
                        structOfClusters(clusterNumber) = proposedCluster;
                        nHyperAcceptances(j) = nHyperAcceptances(j) + 1;
                    end
                    
                end
                
            end
            if(verbose)
                disp('Sampling hyperparameters...')
                disp(['Hypers acceptance rate = ', num2str(nHyperAcceptances./nHyperProposals)]);
            end
        end
    end

    if(mod(sampleNumber,10) == 0 || sampleNumber == nMcmc)
        clear savedDrawFigures savedHyperParameterSamplingFrequency savedNSamples savedVerbose;
        save([saveFileName '.mat']);
    end
    if(verbose)
        disp(['alpha = ', num2str(alpha)]);
        disp(['nClusters =', num2str(length(unique(clusterIDs)))])
%         toc
    end
    
    if(drawFigures)
        refresh
        nClustersOld = doPlots(clusterIDs, data,featureNames, nClustersOld, timeCourseSwitch, multinomialSwitch, bagOfWordsSwitch);
        pause(0.1)
    end

    
end

end


function nClustersOld = doPlots(clusterIDs, data,featureNames, nClustersOld, timeCourseSwitch, multinomialSwitch, bagOfWordsSwitch)

uniqueClusters = unique(clusterIDs);
nClusters = length(uniqueClusters);


if(timeCourseSwitch)
    times = featureNames;
    myCols = cool(nClusters);
    nPlots = ceil(sqrt(nClusters));
    if(nClusters ~= nClustersOld)
        clf;
        nClustersOld = nClusters;
    end
    counter = 1;
    for i = uniqueClusters
        mysubplot(nPlots,nPlots,counter)
        plot(times, data(clusterIDs == i,:), 'color', myCols(counter,:));
        set(gca, 'Color', [.3 .3 .3],'YTickLabel',[],'XTickLabel',[],...
            'YTick',[],'XTick',[], 'xminortick','off', 'yminortick',...
            'off', 'xminorgrid','off', 'yminorgrid','off')
        
        counter = counter + 1;
    end
end
if(multinomialSwitch)
    newData = [];
    for i = uniqueClusters
        newData = [newData; data(clusterIDs == i,:)];
    end
    imagesc(newData');
    if(length(unique(newData(:))) == 3)
        mycolormap = [0 0 1; 1 1 1; 1 0 0];
        colormap(mycolormap);
    else
        colormap jet;
    end
    tickPoints = cumsum(histc(clusterIDs, uniqueClusters));
    
    set(gca,'XTickLabel',[],'YTickLabel',[],'YTick',[],'XTick',[])
    set(gca,'XTick',tickPoints+0.5, 'TickLength', [.4 .4],'LineWidth', 1.5)
end
if(bagOfWordsSwitch)
    
    mygray = gray;
    mygray = mygray(end:-1:1,:);
    
    
    
    
    counter = 1;
    allClusterData = [];
    newData = [];
    for i = uniqueClusters
        tempData = data(clusterIDs == i,:);
        if(mod(counter,2))
            tempData(tempData == 0) = 0.15;
        end
        newData  = [newData; tempData];
        counter = counter + 1;
    end
    
    imagesc(newData')
    colormap(mygray);
    
    tickPoints = cumsum(histc(clusterIDs, uniqueClusters));
    
    set(gca,'XTickLabel',[],'YTickLabel',[],'YTick',[],'XTick',[])
    set(gca,'XTick',tickPoints+0.5, 'TickLength', [.4 .4],'LineWidth', 1)
    
    
    
end
end
