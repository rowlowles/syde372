function imageClassification(data, testData, name, chartTitle)

image_data = extractModels(data);
test_data = extractModels(testData);
[X,Y, classifier] = GEDFilter(image_data);

data = cell(20);
for i = 1:10
   data{(i-1)*2+1} = test_data{(i-1)*4+1};
   data{(i-1)*2+2} = i;
end

classifier = classifyPoints(X, Y, classifier, data);

conf_classifier = confusionmat(classifier(:,1), classifier(:,2));

conf_figure = figure;
confusionchart(conf_classifier);
title(chartTitle);
saveas(conf_figure, "conf_" + name,'png');

end

