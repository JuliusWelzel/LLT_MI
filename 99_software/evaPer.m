% Function that evaluates a binary classifiers performance.
%  
% Input:
%  h: a nx1 binary or probabilistic vector specifying the hypothesized class outcome.
%  y: a nx1 binary vector specifying the actual class.
%  db: an integer specifying the decison boundary (optional)
% 
% Output:
%  per.acc:    The classifiers accuracy (TP+TN)/(TP+TN+FP+FN)
%  per.pre:    The classifiers precision (Positive predictive value): TP/(TP + FP) 
%              --> a low precision indicates many FP
%  per.rec:    The classifiers recall (Sensitivity /True positive rate): TP/(TP+FN) 
%              --> a low recall indicates many FN
%  per.spe:    The classifiers specifity (TNR): TN/(TN+FP)
%  per.f1:     The classifiers F1 value: 2*(precision*recall)/(precision+recall) 
%              --> F1 conveys the balance between precision and recall.
%  per.conf:   A confusion matrix (
%              ---> columns: actual class labels, rows: predicted class labels
%  per.trials: An integer documenting the number of trials
%  per.sig:    A boolean documenting whether tge accuracy was significant.
%
%  per.brier:  The classifiers brier accuracy (if probabilistic values are provided)
%
% Example usage:
%  h = [1 1 1 1 1 0 0 0 0 0];
%  y = [1 1 1 0 0 0 0 0 0 1];
%  eva = evaPer(h,y)
%
% Niclas Braun, 21.10.2015

function per = evaPer(h,y,db)

% Check whether default decision boundary shall be used
if nargin <3
    db = 0.5;
end

% Transform probabilistic values into binary values
hb = h;
hb(h>=db) = 1; hb(h<db) = 0;

%% Calculate confusion matrix:
P = hb(find(y==1)); 
N = hb(find(y==0));

TP = sum(P); per.con(1,1)  = TP;
FP = sum(N); per.con(1,2)  = FP;
FN = sum(~P); per.con(2,1) = FN;
TN = sum(~N); per.con(2,2) = TN;

%% Calculate performance measures
per.acc = (TP+TN)/(TP+TN+FP+FN);
per.pre = TP/(TP + FP); 
per.rec= TP/(TP+FN);
per.spe= TN/(TN+FP);
per.f1=  2*(per.pre * per.rec)/(per.pre+per.rec);
per.trials = length(y);

%% Check for significance, if check_classifier.m is available
if exist('check_classifier') == 2
    [X,chanceLevel] = check_classifier([length(y==0) length(y==1)], 0.05, 2);
    if per.acc > chanceLevel
       per.sig = 1;
    else
        per.sig = 0;
    end
end

%% Calculate brier, if probabilistic values are provided
if (sum(h==0) && sum(h==1)) ~= length(h)
    per.brier = (1/length(h)) * sum((h-y).^2);
end
    

