function setup_clusters( wokersnum )
% Set up the number of wokrers
% if wokersnum is greater than maximum of cores,
% set the number of wokrers to be maximum of cores-1

% Get the number of workers
% myCluster = parcluster('local');
%  display(myCluster.NumWorkers);

% Setup for parallel computing

% Get the number of cores
numcores = feature('numcores');

if numcores > wokersnum
    parpool('local',wokersnum);
else
    parpool('local',max(numcores-1,1));
%     parpool('local',2);
end

% Shut down the parallel pool
% delete(gcp('nocreate'));