function [X,Y, classifier]  = nearestNeighbourFilter(k, varargin)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

% calculate minValue
numClasses = size(varargin,2)/ 2;
mins = zeros(numClasses, 2);
maxs = zeros(numClasses, 2);
for i = 1:numClasses
    class = varargin{(i-1)*2+1};
    mins(i,:) = min(class);
    maxs(i,:) = max(class);
end

minValue = floor(min(mins))-5;
maxValue = ceil(max(maxs))+5;

feature1Vals = minValue(1):0.05:maxValue(1);
feature2Vals = minValue(2):0.05:maxValue(2);

[X, Y] = meshgrid(feature1Vals, feature2Vals);

arrSize = [size(feature2Vals,2) size(feature1Vals,2)];

classDist = zeros(numClasses,arrSize(1),arrSize(2));

for i = 1:size(feature1Vals,2)
    for j = 1:size(feature2Vals,2)
        for class = 1:numClasses
            classData = varargin{(class-1)*2+1};
            [~, index] = findNearestNeighbours(classData, [feature1Vals(1,i) feature2Vals(1,j)], k);
            prototype = mean(classData(index, :),1);
            dist = sqrt(sum(([feature1Vals(1,i) feature2Vals(1,j)]-prototype).^2));
            classDist(class,j,i) = dist;
        end
    end
end

classifier = zeros(arrSize);
for class = 1:numClasses
   classifier(classDist(class,:,:) < min(classDist(1:end ~= class,:,:),[],1)) = class;
end

figure;
plotClasses(varargin{:});
hold on;
contour(X,Y,classifier);

end

