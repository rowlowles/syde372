% Lab 1
close all;
clear all;
% Generate clusters
% Class A
nA = 200;
muA = [5;10];
sigmaA = [8 0; 0 4];
classA = repmat(muA',[nA, 1])  + randn(nA,2)*chol(sigmaA);
meanA = mean(classA);

% Class B
nB = 200;
muB = [10;15];
sigmaB = [8 0; 0 4];
classB = repmat(muB',[nB, 1]) + randn(nB,2)*chol(sigmaB);
meanB = mean(classB);

% Class C
nC = 100;
muC = [5;10];
sigmaC = [8 4; 4 40];
classC = repmat(muC',[nC, 1]) + randn(nC,2)*chol(sigmaC);
meanC = mean(classC);


% Class D
nD = 200;
muD = [15;10];
sigmaD = [8 0; 0 8];
classD = repmat(muD',[nD, 1]) + randn(nD,2)*chol(sigmaD);
meanD = mean(classD);

% Class E
nE = 150;
muE = [10;5];
sigmaE = [10 -5; -5 20];
classE = repmat(muE',[nE, 1]) + randn(nE,2)*chol(sigmaE);
meanE = mean(classE);


% Plot the distributions' means
% Plot distribution A & B and STD contour
figure;
plotClasses(classA,'Class A',classB,'Class B');
hold on; 
plotStdContours([1], meanA, sigmaA, meanB, sigmaB);

% Plot distribution C, D & E and STD contour
figure;
plotClasses(classC,'Class C',classD,'Class D', classE,'Class E');
hold on; 
plotStdContours([1], meanC, sigmaC, meanD, sigmaD, meanE, sigmaE);




