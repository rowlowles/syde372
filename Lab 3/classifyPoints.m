function [elements] = classifyPoints(X,Y, classifier, varargin)
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here
    varargin = varargin{1};
    classNum =  size(varargin,2)/ 2;
    elements = [];
    
    for i = 1:classNum
       originalPoints = varargin{(i-1)*2+1};
       for index = 1:size(originalPoints, 1)
          point = originalPoints(index,:);
          point;
          classification = findClassification(X,Y,classifier,point(1,1), point(1,2));
          elements = [elements; i classification];
       end
    end
end

