function C = uniquePairs(A,B)
% Finds how many times a pair of values occurs

%A = [1 6 1 10 1 7 1 9  3 6 0 0 -2 -2];
%B = [7 2 3  5 6 8 7 9 10 2 0 0 1 0];

%// Find the unique pairs
AB  = [A;B].';
ABu = unique(AB, 'rows');

%// Count the number of occurrences of each unique pair
%// (pretty sure there's a better way to do this...)     
C = arrayfun(@(ii) ...
    [ABu(ii,:) sum(all(bsxfun(@eq, AB, ABu(ii,:)),2))], 1:size(ABu,1), ...
    'UniformOutput', false);
C = cat(1,C{:});