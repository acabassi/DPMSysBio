function doDiagnosticPlots(fileName, uniqueIdentifier)
saveFileName = [strtok(fileName, '.'),'_Results_Chain', num2str(uniqueIdentifier)];
    
load(saveFileName)

allData = csvread([saveFileName '.csv'], 1,0);

alphas   = allData(:,1);
clusters = allData(:,2:(nGenes+1));

if(strcmp(dataType, 'TimeCourseUnequalSpacing'))
   hypers = allData(:,(nGenes+1):end);  
end

sorted = sort(clusters,2).' ;
temp   = [true(1,size(sorted,2)) ; diff(sorted)>0] ;
nClusters = sum(temp);

figure
subplot(2,2,1)
plot(alphas)
xlabel('Number of iterations')
ylabel('alpha (concentration parameter)')
subplot(2,2,2)
plot(nClusters)
xlabel('Number of iterations')
ylabel('Number of clusters')
subplot(2,2,3)
hist(alphas)
xlabel('alpha (concentration parameter)')
subplot(2,2,4)
hist(nClusters)
xlabel('Number of clusters')

end