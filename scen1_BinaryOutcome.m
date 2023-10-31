function results = scen1_BinaryOutcome(seed,n,run_ID,isdisplay)
%% The source code is revised from the R code from
%  Roy, J., Lum, K.J., Zeldow, B., Dworkin, J.D., Re III, V.L.
% and Daniels, M.J., 2017. Bayesian nonparametric generative models
% for causal inference with missing at random covariates. Biometrics.

% tic

%% Set random seed-------------------------------------------------------------
rng(seed);

inputfile = ['datasets/scen1/scen1_n'  num2str(n)   '_run_ID'  num2str(run_ID) '.txt'];
tmpdata =  dlmread(inputfile);
pt = 1;
p1 = 2;
p2 = 2;
T = tmpdata(:,pt);
X = tmpdata(:,pt+1:end-1);
Y = tmpdata(:,end);

% Standardize continuous covariates (important when choosing priors).
X(:,p1+1:p1+p2) = zscore(X(:,p1+1:p1+p2));

% Define design matrix.
X_matrix = [ones(size(T,1),1),X];
TX = [T,X];
TX_matrix = [ones(size(T,1),1),TX];

% End  generation of dataset-----------------------------------------------

%% Define constants--------------------------------------------------------------
if n < 1000 
    % Define number of Markov Chain Monte Carlo (MCMC) draws in Gibbs Sampler.

    % For small n (e.g. n = 250), need more samples.
     GIBBS_TOTAL = 10000;

    % Define number of MCMC draws to 'burn' in Gibbs Sampler (check convergence).

    % For small n (e.g. n = 250).
      BURN_IN = 2000; 

else
    % Define number of Markov Chain Monte Carlo (MCMC) draws in Gibbs Sampler.
    
    % For big n (e.g. n = 1000), can use less samples (check posterior).
    GIBBS_TOTAL = 10000;

    % Define number of MCMC draws to 'burn' in Gibbs Sampler (check convergence).
    BURN_IN = 2000;
end

%  The value of auxiliary parameters for new clusters in Gibbs
%  sampling.
NUM_AUX_PARAM = 5;

% Define number of Monte Carlo draws when using Monte Carlo integration to
% integrate out confounders.
NUM_MC = 10000;

% Define number of observations per draw when using Monte Carlo integration to
% integrate out confounders.
M = 1000;

% End Define constants --------------------------------------------------


%% Set initial values ------------------------------------------------------------

%% Hyperparameters-----------------------------------------------------------------

% For modeling binary covariates.
a0 = 1;
b0 = 1;

% For modeling continuous covariates.
nu0 = 2;
tau0 = 1;
c0 = 0.5;
mu0 = 0;

% For concentration parameters.
alpha_shape = 1;
alpha_rate = 1;

%% Initialize alpha parameters-----------------------------------------------------
alpha_theta = 1;
alpha_psi = 1;
alpha_omega = 1;

%% Parameters for cluster memberships------------------------------------------

% Set initial values for cluster membership,y,t,x
% Use k-means to find  data clusters
z = init_cluster_kmeans(Y,TX,X,2,2,2);

% % Calculate number of clusters
[k,kjs,kljs] = computeclusnum(z);

% End Set initial values.

%% CHANGE OCCURS HERE------------------------------------------------------------

%% Make file_name where results will be stored for each dataset-------------------
resultsfile = ['results/scen1/cedp_scen1_n'  num2str(n)   '_run_ID'  num2str(run_ID) '.mat'];
workplacefile = ['results/scen1_workplace/cedp_scen1_n'  num2str(n)   '_run_ID'  num2str(run_ID) '.mat'];

%% Update values of parameters-----------------------------------------------------
% % For outcome regression coefficients parameters in outcome regression.
[beta_coeff0,diagbetacov0,beta_coeff] = init_beta_params_bern(Y,TX,TX_matrix,z,k);

%% Update values of parameters for treatment regression coefficients-----------------
[gamma_coeff0,diaggammacov0,gamma_coeff] = init_gamma_params(T,X,X_matrix,z,k,kjs);

%% Update covariate parameters------------------------------------------------------
% Within each cluster parameters are stored in a long vector nXClus x 1.
% Binary covariates have 1 parameter (x_pi_param).
% Continuous covariates have 2 parameters (x_mean_param, x_var_param).
x_pi_param = update_binarycovs_params(X,p1,z,k,kjs,kljs,a0,b0);

[x_mean_param,x_var_param] = update_continuouscovs_params(X,p1,p2,z,k,kjs,kljs,nu0,tau0,c0,mu0);

%% For MC integration------------------------------------------------------------
% Calculate f0y,f0t,f0x for use in calculating causal effect in Gibbs Sampler.
% Average covariate distribution over prior for each x_i
% variables of x.

beta_coeff_prior_mc = mvnrnd(beta_coeff0,diagbetacov0,NUM_MC);

gamma_coeff_prior_mc = mvnrnd(gamma_coeff0,diaggammacov0,NUM_MC);

% Make vectors to store draws from Gibbs Sampler.
psi_rr = NaN(GIBBS_TOTAL,1);
psi_rd = NaN(GIBBS_TOTAL,1);

%% Parameters for analysing----------------------------------------------------------

% Store the number of each y cluster.
yclusNum = zeros(GIBBS_TOTAL,1);

% Store the  log likehood of data.
loglikeregData = zeros(GIBBS_TOTAL,1);

% Store the values of alpha.
alpha_theta_draws = zeros(GIBBS_TOTAL,1);
alpha_psi_draws = zeros(GIBBS_TOTAL,1);
alpha_omega_draws = zeros(GIBBS_TOTAL,1);

%% Start Gibbs Sampler -------------------------------------------------------
for gibbsreps = 1:GIBBS_TOTAL
    % First draw each parameter for BNP model. Then calculate causal effect.
    
    %% Update cluster membership-------------------------------------------------------
    [z,beta_coeff,gamma_coeff,x_pi_param,x_mean_param,...
        x_var_param] = update_cluster_binary(T,X,Y,p1,p2,z,beta_coeff,gamma_coeff,x_pi_param,x_mean_param,...
        x_var_param,alpha_theta,alpha_psi,alpha_omega,NUM_AUX_PARAM,beta_coeff0, diagbetacov0,...
        gamma_coeff0,diaggammacov0,c0, mu0, nu0, tau0, a0, b0);
    
   if mod(gibbsreps,1000) == 0
      fprintf([' sample ' num2str(gibbsreps),' k ', num2str( size(beta_coeff,1)) '\n']);
   end
    % End Update cluster membership--------------------------------------------------
    
    %% Calculation of parameters of clusters------------------------------------
    %  Calculate number of clusters.
    [k,kjs,kljs] = computeclusnum(z);
    
    % Calculate numbers of data in each cluster.
    [njs,nljs,nhljs] = computenuminclus(z,k,kjs,kljs);
    
    %% Update parameters------------------------------------------------------------------------
    % Update beta parameters.
    beta_coeff = update_betas_bern(Y,TX_matrix,z,k,beta_coeff, beta_coeff0,diagbetacov0);
    
    % Update parameters for T.
    gamma_coeff = update_gammas(T,X_matrix,z,k,kjs,gamma_coeff,gamma_coeff0,diaggammacov0);
    
    % Update parameters for X.
    if p1 > 0
        x_pi_param = update_binarycovs_params(X,p1,z,k,kjs,kljs,a0,b0);
    end
    
    % Update parameters for continuous covariates.
    if p2 > 0
        [x_mean_param,x_var_param] = update_continuouscovs_params(X,p1,p2,z,k,kjs,kljs,nu0,tau0,c0,mu0);
    end
    
    % End Update parameters --------------------------------------------
    
    %% Store the number of each y cluster  and calculate the log likehood of data---------------
    yclusNum(gibbsreps) = k;
    loglikeregData(gibbsreps) = computeloglike_binary(T,X,Y,TX_matrix,X_matrix,p1,p2,...
        z,k,kjs,kljs,beta_coeff,gamma_coeff,x_pi_param,x_mean_param,x_var_param);
    
    %% Update concentration parameters--------------------------------------------
    alpha_theta =  update_alpha_theta(alpha_shape,alpha_rate,n,k,alpha_theta);
    alpha_psi = update_alpha_psi(alpha_shape,alpha_rate,z,alpha_psi);
    alpha_omega =  update_alpha_omega(alpha_shape,alpha_rate,z,alpha_omega);
    
    alpha_theta_draws(gibbsreps,1) = alpha_theta;
    alpha_psi_draws(gibbsreps,1) = alpha_psi;
    alpha_omega_draws(gibbsreps,1) = alpha_omega;
    
    % End Update concentration parameters-----------------------------------
    
    % End Update of all parameters in BNP model--------------------------------
    
    %% Calculate causal effect using Monte Carlo Integration to integrate out confounders.
    % Don't need to do every iteration so compute every 100th draw after
    % BURN_IN).
    % Start----------------------------------------------------------------------
    
     if gibbsreps >= BURN_IN && mod(gibbsreps,100) == 0
        %% Calculation of causal effects
        [psi_rr(gibbsreps,1),psi_rd(gibbsreps,1)] = compute_effect_binary(X,n,p1,p2,z,k,kjs,kljs,njs,nljs,nhljs,...
            beta_coeff,gamma_coeff,x_pi_param,x_mean_param,x_var_param,...
            alpha_theta,alpha_psi,alpha_omega,beta_coeff0,diagbetacov0,gamma_coeff0,diaggammacov0,...
            c0, mu0, nu0, tau0, a0, b0,M,beta_coeff_prior_mc,gamma_coeff_prior_mc);
        
        % Error handling
        if isfinite(psi_rr(gibbsreps,1)) == false
            return;
        end
        
    end
    
    
end  % End Gibbs Sampler loop ---------------------------------------------

%% Calculate summary statistics from posterior distributions
median_alpha_theta = median(alpha_theta_draws(BURN_IN:GIBBS_TOTAL,1));
median_alpha_psi = median(alpha_psi_draws(BURN_IN:GIBBS_TOTAL,1));
median_alpha_omega = median(alpha_omega_draws(BURN_IN:GIBBS_TOTAL,1));


% True values of average relative risk (for this scenario)
psi_rr_true = 1.5;

[median_psi_rr,sd_psi_rr,cov_psi_rr,width_psi_rr] = compute_psi_rd_stats(psi_rr_true,psi_rr);


% True values of average risk difference (for this scenario)
psi_rd_true = 0.13;

[median_psi_rd,sd_psi_rd,cov_psi_rd,width_psi_rd] = compute_psi_rd_stats(psi_rd_true,psi_rd);

% Save results to a vector
results = [median_alpha_theta,median_alpha_psi,median_alpha_omega,median_psi_rr,sd_psi_rr,cov_psi_rr,width_psi_rr,median_psi_rd,sd_psi_rd,cov_psi_rd,width_psi_rd];

% Output results to a file
save(resultsfile,'results');
save(workplacefile);

%% Visualize trace plot for analysing

if isdisplay == 1
% Plot y clusters
plot_yclusNum_hist(yclusNum);
plot_yclusNum(yclusNum);

% Plot log likelihood of data
plot_loglike(loglikeregData);

% Plot concentration parameters
plot_alpha_theta(alpha_theta_draws);
plot_alpha_psi(alpha_psi_draws);
plot_alpha_omega(alpha_omega_draws);

% Plot relative causal risk and average causal effects
plot_psi_rr(psi_rr(~isnan(psi_rr)));
plot_psi_rd(psi_rd(~isnan(psi_rd)));

end

% toc