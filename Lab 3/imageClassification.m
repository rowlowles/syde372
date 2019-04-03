function [classifier, X, Y] = imageClassification(data, testData, name, chartTitle)

image_data = extractModels(data);
test_data = extractModels(testData);
[X,Y, classifier] = GEDFilter(image_data);

data_cell = cell(20);
for i = 1:10
   data_cell{(i-1)*2+1} = test_data{(i-1)*4+1};
   data_cell{(i-1)*2+2} = i;
end

classification = classifyPoints(X, Y, classifier, data_cell);

conf_classifier = confusionmat(classification(:,1), classification(:,2));

conf_figure = figure;
confusionchart(conf_classifier);
title(chartTitle);
saveas(conf_figure, "conf_" + name,'png');

end

