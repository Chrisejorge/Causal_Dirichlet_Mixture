function [expitTable,MAX_EXP,EXP_TABLE_SIZE] = generate_expitTable()

MAX_EXP = 6; % [-MAX_EXP,+MAX_EXP], if the value x < -MAX_EXP, take  expit(x)
% as expit(-MAX_EXP), x > MAX_EXP, take  expit(x) as expit(MAX_EXP)
EXP_TABLE_SIZE = 1000; % cut [-MAX_EXP,MAX_EXP] into EXP_TABLE_SIZE intervals

expitTable = 1./(1 + exp(-(-MAX_EXP + (0:EXP_TABLE_SIZE-1)*(2 * MAX_EXP)/EXP_TABLE_SIZE)));
