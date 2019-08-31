% find closest value and its idx of a specified in an array
function [idx, value] = findC(val,array)
tmp = abs(array-val);
[tmpMin idx] = min(tmp); %index of closest value
value = array(idx); %closest value



