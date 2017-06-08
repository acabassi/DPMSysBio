fHandle          = @DPCluster_TimeCourseUnequal;  % Likelihood function to use for time series
fileName         = 'data-30-2.csv';       % Data file
uniqueIdentifier = 1;    % ID for output file (change if you re-run, to avoid overwriting old results)
nSamples         = 20000;  % Number of MCMC iteration
dataType         = 'TimeCourseUnequalSpacing';  % Leave this!
drawFigures      = false;  % Should figures be plotted ?
hyperParameterSamplingFrequency = 1;  % Leave this (unless very slow, in which case try increasing to 2, then 5, then 10)
verbose          = false; % Whether or not to print output to screen
initialise       = true;  % If you have to stop the sampler, and then wish to rerun at a later date, set this to false 
gammaPrior       = [2 4]; % Prior for the GP hyperparameter (leave this)
thinningFreq     = 10;     % If set to X, records every X-th sample
inputSeed        = 1;     % Set to NaN if you do not want a seed
GUIflag          = false; % Leave this

feval(fHandle, fileName, uniqueIdentifier, nSamples, dataType, ...
            drawFigures, hyperParameterSamplingFrequency, verbose, initialise, ...
            gammaPrior, thinningFreq,inputSeed,GUIflag);
