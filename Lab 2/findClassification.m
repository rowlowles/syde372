function [classNum] = findClassification(X,Y,classifier, x, y)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here

X_inc = abs(X(1,1) - X(end, end))/(length(X)-1);
Y_inc = abs(Y(1,1) - Y(end, end))/(length(Y)-1);

X_indx = floor((x-X(1,1))/X_inc)+1;
Y_indx = floor((y-Y(1,1))/Y_inc)+1;

% Safe guard for out of range values
if Y_indx <= 0 
    Y_indx = 1;
end

if X_indx <= 0 
    X_indx = 1;
end

if Y_indx > size(classifier,1) 
    Y_indx = 1;
end

if X_indx > size(classifier,2) 
    X_indx = 1;
end


classNum = classifier(Y_indx, X_indx); 

end

