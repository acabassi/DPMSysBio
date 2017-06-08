clear
close all
randn('seed', 100)
rand('seed', 100)

%%% Run DPMSysBio for time course data       
fileName  = 'TimeCourseExample.csv';

samplingFreq     = 1;
nSamples         = 10;
dataType         = 'TimeCourse';
verbose          = true;
initialise       = false;
drawFigures      = true;
uniqueIdentifier = 1;
inputSeed        = 100;
gammaPrior       = [2 4];
hyperParameterSamplingFrequency = 1;

DPM(fileName, uniqueIdentifier, nSamples, dataType, ...
    drawFigures, hyperParameterSamplingFrequency, verbose, initialise,...
    gammaPrior, samplingFreq, inputSeed)