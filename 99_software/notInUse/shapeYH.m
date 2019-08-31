% auxiliary function that turns the probabilistic classifier outcomes
% into a binary format (1s for class 1, 2s for class2)
% Probabilistic outcomes are interpreted relative to class 1
% (e.g. [0.2 0.7 0.7]' --> [2 1 1])

function yh = shapeYH(yh)

for i=1:length(yh)
   if yh(i) > 0.5
            yh(i) = 1;
   else
            yh(i) = 2;
   end      
end


