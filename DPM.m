function DPM(fileName, uniqueIdentifier, nSamples, dataType, ...
    drawFigures, hyperParameterSamplingFrequency, verbose, initialise,...
    gammaPrior, thinningFreq, inputSeed)


switch dataType
    case 'TimeCourse'
        fHandle = @DPCluster_TimeCourseUnequal;
        dataType = 'TimeCourseUnequalSpacing';    
    case 'TimeCourseUnequal'
        fHandle = @DPCluster_TimeCourseUnequal;
        dataType = 'TimeCourseUnequalSpacing';    
    case 'TimeCourseEqual'
        fHandle = @DPCluster_TimeCourseUnequal;
        dataType = 'TimeCourseUnequalSpacing';
    case 'GPR'
        fHandle = @DPCluster_TimeCourseUnequal;
        dataType = 'TimeCourseUnequalSpacing';
    case 'Multinomial'
        fHandle = @DPCluster_Multinomial;
        dataType = 'Multinomial';
    case 'DiscreteGeneExpression'
        fHandle = @DPCluster_Multinomial;
        dataType = 'Multinomial';
    case 'Discrete'
        fHandle = @DPCluster_Multinomial;
        dataType = 'Multinomial';
    case 'BagOfWords'
        fHandle = @DPCluster_BagOfWords;
        dataType = 'BagOfWords';
    case 'Binary'
        fHandle = @DPCluster_BagOfWords;
        dataType = 'BagOfWords';
        
end

feval(fHandle, fileName, uniqueIdentifier, nSamples, dataType, ...
    drawFigures, hyperParameterSamplingFrequency, verbose, initialise, ...
    gammaPrior, thinningFreq,inputSeed, false);


end