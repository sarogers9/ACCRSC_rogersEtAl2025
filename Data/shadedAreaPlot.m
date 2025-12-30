function [mu,sem] = shadedAreaPlot(data,color)

%this function makes bar graphs with SEM error bars
%INPUTS =   data - your data matrix, with variables to plot in columns and
%                  observations/replicates in rows
%OUTPUTS =  mu - mean
%           sem - standard error
%           bar plot
%
%Written by Sophie A. Rogers, Corder Laboratory, University of Pennsylvania
%%
mu = nanmean(data,2);
sem = nanstd(data,[],2)./sqrt(size(data,2)-1);
plot(mu,'color',color)
hold on
y_upper = mu+sem;
y_lower = mu-sem;
fill([1:size(data, 1), fliplr(1:size(data, 1))], [y_upper; flipud(y_lower)], color, 'FaceAlpha', 0.2, 'EdgeColor', 'none');