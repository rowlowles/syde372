function [X,Y, classifier]  = GEDFilter(varargin)
%GEDFilter Summary of this function goes here
%   varargin: class information (any number of classes >=2 can be included)
%             1st class arg: class data
%             2nd class arg: class mean
%             3rd class arg: class covariance matrix
%             4th class arg: class name

% calculate minValue
varargin = varargin{1};
numFieldsPerClass = 4;
numClasses = size(varargin, 1) / numFieldsPerClass;
mins = zeros(numClasses, 2);
maxs = zeros(numClasses, 2);
for i = 1:numClasses
    class = varargin{(i-1)*4 + 1};
    mins(i,:) = min(class);
    maxs(i,:) = max(class);
end

minValue = min(mins);
maxValue = max(maxs);

iteration_amount = (maxValue-minValue)/750;
minValue = minValue - 2*iteration_amount;
maxValue = maxValue + 2*iteration_amount;

feature1Vals = minValue(1):iteration_amount(1):maxValue(1); % step was 0.05
feature2Vals = minValue(2):iteration_amount(2):maxValue(2);

[X, Y] = meshgrid(feature1Vals, feature2Vals);

arrSize = [size(feature2Vals,2) size(feature1Vals,2)];

classDist = zeros(numClasses,arrSize(1),arrSize(2));

% Iterate through distance grid
for i = 1:size(feature1Vals,2)
    for j = 1:size(feature2Vals,2)
        for class = 1:numClasses
            % Extract class info from varargs
            classMean = varargin {(class-1)*numFieldsPerClass + 2};
            classCov = varargin {(class-1)*numFieldsPerClass + 3};
            
            point = [feature1Vals(1,i) ;feature2Vals(1,j)];
            classDist(class,j,i) = sqrt(sum((point-classMean).'*inv(classCov)*(point-classMean)));
        end
    end
end

classifier = zeros(arrSize);
for class = 1:numClasses
   classifier(classDist(class,:,:) < min(classDist(1:end ~= class,:,:),[],1)) = class;
end

% Plot the data
figure
classDataAndNameIdxs = zeros(1, numClasses*2);
idx = 1;
for classIdx = 1:numFieldsPerClass:length(varargin)
    classDataAndNameIdxs(idx) = classIdx;
    classDataAndNameIdxs(idx+1) = classIdx + 3;
    idx = idx + 2;
end
plotArgs=varargin(classDataAndNameIdxs);
plotClasses(plotArgs{:});
hold on;
contour(X,Y,classifier, 'DisplayName', 'GED Boundary');

end
