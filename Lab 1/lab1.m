% Lab 1
close all;
clear all;
% Generate clusters
% Class A
nA = 200;
muA = [5;10];
sigmaA = [8 0; 0 4];
classA = repmat(muA,nA/2,2) + randn(nA,2)*chol(sigmaA);
meanA = mean(classA,2);

% Class B
nB = 200;
muB = [10;15];
sigmaB = [8 0; 0 4];
classB = repmat(muB,nB/2,2) + randn(nB,2)*chol(sigmaB);
meanB = mean(classB,2);

% Class C
nC = 100;
muC = [5;10];
sigmaC = [8 4; 4 40];
classC = repmat(muC,nC/2,2) + randn(nC,2)*chol(sigmaC);
meanC = mean(classC,2);


% Class D
nD = 200;
muD = [15;10];
sigmaD = [8 0; 0 8];
classD = repmat(muD,nD/2,2) + randn(nD,2)*chol(sigmaD);
meanD = mean(classD,2);

% Class E
nE = 150;
muE = [10;5];
sigmaE = [10 -5; -5 20];
classE = repmat(muE,nE/2,2) + randn(nE,2)*chol(sigmaE);
meanE = mean(classE,2);

% Plot the distributions' means
% TODO: What is a unit standard deviation contour? Ask during tutorial or
% see when it comes up in lecture. 


