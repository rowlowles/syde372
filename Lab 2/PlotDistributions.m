function PlotDistributions(estimatedDistribution,actualDistribution, rangeOfX, chartTitle, fileName)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
Y_a_d = actualDistribution(rangeOfX);
Y_e_d = estimatedDistribution(rangeOfX);

currFig = figure;
plot(rangeOfX, Y_a_d)
hold on
plot(rangeOfX, Y_e_d)
legend({'actual','estimate'})
title(chartTitle)
xlabel('Value') 
ylabel('Probability') 

saveas(currFig,fileName,'png')
end

