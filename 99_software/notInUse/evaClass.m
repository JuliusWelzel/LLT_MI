function res = evaClass(y,yh)
% Function that evaluates the performance of a binary classifier. 
% Provided are the error rate, sensitivity, specifity, positive predictive value 
% and negative predictive value.
%
% Input:
% y = a nx1 vector of target values (1s for 1st class, 2s for 2nd class)
% yh = a nx1 vector of actual outcomes (1s for 1st class, 2s for 2nd class)
%
% Output:
%
% Example usage:  
% y = [1 1 1 2 2 2]'; yh= [1 1 1 2 2 1]'; res = contingency_table(y,yh)
% (adapted from: % http://matlabgeeks.com/tips-tutorials/clustering-part-1-binary-classification/)
%keyboard;

% True positives
TP = yh(y == 2); % find all targets of class 2
TP = length(TP(TP == 2)); % count those that have been classified as class 2

% False positives
FP = yh(y == 1); % find all targets of class 1
FP = length(FP(FP == 2));  % count those that have been classified as class 2 

% False negatives
FN = yh(y == 2); % find all targets of class 2
FN = length(FN(FN == 1)); % count those that have been classified as class one

% True negatives
TN = yh(y == 1); % find all targets of class 1
TN = length(TN(TN == 1)); % count all that have been classified as class 1

% Sensitivity
res.sen = TP/(TP+FN);

% Specicfity
res.spec = TN/(FP+TN);

% Accuracy
%keyboard;
res.acc = (TP+TN) / (TP+FN+FP+TN); 