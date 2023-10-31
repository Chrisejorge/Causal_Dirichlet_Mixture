function logitinv =  expit2(x)
% Define function to calculate inverse of the logit function
% x column vector

MAX_EXP = 6; % [-MAX_EXP,+MAX_EXP], if the value x < -MAX_EXP, take  expit(x)
% as expit(-MAX_EXP), x > MAX_EXP, take  expit(x) as expit(MAX_EXP)
EXP_TABLE_SIZE = 1000; % cut [-MAX_EXP,MAX_EXP] into EXP_TABLE_SIZE intervals

if x > MAX_EXP
    logitinv = expitTable(EXP_TABLE_SIZE);
elseif x < -MAX_EXP
    logitinv = expitTable(1);   
else
    logitinv = expitTable(int32((x + MAX_EXP) * EXP_TABLE_SIZE / (2 * MAX_EXP)  + 1));
end
