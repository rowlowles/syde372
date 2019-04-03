function [cellOutput] = extractModels(data)
% row 1 - feature 1
% row 2 - feature 2
% row 3 - feature 3

%             1st class arg: class data
%             2nd class arg: class mean
%             3rd class arg: class covariance matrix
%             4th class arg: class name
vargIndex = 1;

classes = unique(data(3,:));
cellOutput = cell(length(classes)*4,1);

for class = classes
    extractedClassData = data(:, find(data(3,:) == class));
    cellOutput{vargIndex} = extractedClassData(1:2,:)';
    cellOutput{vargIndex+1} = mean(extractedClassData(1:2,:),2);
    cellOutput{vargIndex+2} = cov(extractedClassData(1:2,:)');
    cellOutput{vargIndex+3} = num2str(class);  
    vargIndex = vargIndex + 4;
end

end

