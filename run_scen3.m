% Init path for matlab
initpath

% Set up clusters for parallel computing
% setup_clusters(5);

%% Set random seed-------------------------------------------------------------
load('seedfile.mat','seed')

% Define number of observations for each dataset.
% n = 250;
n = 1000;

% Variable to decide to visualize trace plot for analysing, 1 yes, 0 no
isdisplay = 0;

% The number of  datasets
NUM_DATA = 10;

tic

Error=zeros(1,NUM_DATA);

parfor run_ID = 1:NUM_DATA
   try
    scen3_results(run_ID,:) = scen3_ContinuousOutcome(seed,n,run_ID,isdisplay);
    Error(run_ID)=0;
  catch ME
    disp(ME.message); %Will show the error info
    Error(run_ID)=1;
  end
end

save(['scen3_n',num2str(n),'_results.mat'],'scen3_results');

toc

% Shut down the parallel pool
delete(gcp('nocreate'));

quit;
