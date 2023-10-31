function logitinv =  expit(x)
% Define function to calculate inverse of the logit function
% x column vector

logitinv = 1./(1 + exp(-x));

% end