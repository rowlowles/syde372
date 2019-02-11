function [dist, index] = findNearestNeighbours(points,searchPoint,k)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
v = points-searchPoint;
dist = sqrt(sum(v.^2,2));

[dist, index] = sort(dist);
dist = dist(1:k);
index = index(1:k);
end

