function [ mean_data sem_data ] = mean_SEM( data )
%mean_SEM calculates the standard error and mean of a vector of data
% Author: Julius Welzel, University of Oldenburg 2018

mean_data = [];
sem_data = [];

if ~iscell(data)
    data = num2cell(data);
end

for r = 1:size(data,1)
    mean_data(r) = mean(data{r},'omitnan');

    sem_data(r) = std(data{r},'omitnan')/sqrt(length(data{r}));
end

end

