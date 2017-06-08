dataType = 'Multinomial';
fileName = 'CellCycleData.csv';

dataType = 'BagOfWords';
fileName = 'HarbisonProcessedUsingOrlandoWildTypeSystematicnames.csv';

dataType = 'TimeCourseUnequalSpacing';
fileName = 'TimeCourseData.csv';


dataType = 'TimeCourseEqualSpacing';
fileName = 'TimeCourseData.csv';


dataType = 'Multinomial';
fileName = 'CellCycleData.csv';


dataType = 'TimeCourseEqualSpacing';
fileName = '27c.csv';


uniqueIdentifier = 1;
nSamples         = 1000;
drawFigures      = true;

hyperParameterSamplingFrequency = 1;
verbose       = true;
initialise    = true;
gammaPrior    = [2 4];
thinningFreq  = 5;

switch dataType
    case 'TimeCourseEqualSpacing'
        fHandle = @DPCluster_TimeCourseEqual;
    case 'TimeCourseUnequalSpacing'
        fHandle = @DPCluster_TimeCourseUnequal;
    case 'Multinomial'
        fHandle = @DPCluster_Multinomial;
    case 'BagOfWords'
        fHandle = @DPCluster_BagOfWords;
end



feval(fHandle, fileName, uniqueIdentifier, nSamples, dataType, ...
    drawFigures, hyperParameterSamplingFrequency, verbose, initialise, ...
    gammaPrior, thinningFreq);
 