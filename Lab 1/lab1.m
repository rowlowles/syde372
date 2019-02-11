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
clusters1 = figure;
plotClasses(classA,'Class A',classB,'Class B');
hold on; 
plotStdContours([1], meanA, sigmaA, meanB, sigmaB);

% Plot distribution C, D & E and STD contour
clusters2 = figure;
plotClasses(classC,'Class C',classD,'Class D', classE,'Class E');
hold on; 
plotStdContours([1], meanC, sigmaC, meanD, sigmaD, meanE, sigmaE);


%% MED Classifier
% Step 1: Find distance between two points
% 
% figure(clusters1)
% % g(x) = [a b c] in form aX_1 + bX_2 + c = 0
% gx = [(muA - muB)', .5*(muB'*muB - muA'*muA)];
% vec = [-gx(2)/gx(1) -gx(3)/gx(1)];
% refline(vec(1), vec(2))
% 
% figure(clusters2)
% gx = [(muC - muE)', .5*(muE'*muE - muC'*muC)];
% vec = [-gx(2)/gx(1) -gx(3)/gx(1)];
% refline(vec(1), vec(2))
% 
% gx = [(muD - muE)', .5*(muE'*muE - muD'*muD)];
% vec = [-gx(2)/gx(1) -gx(3)/gx(1)];
% refline(vec(1), vec(2))
%% MED For Clusters 1
minValue = floor(min(min(classA, classB)));
maxValue = ceil(max(max(classA, classB)));

feature1Vals = minValue(1):0.1:maxValue(1);
feature2Vals = minValue(2):0.1:maxValue(2);

arrSize = [size(feature1Vals,2) size(feature2Vals,2)];
boundary1 = [];

for x1MED = 1:arrSize(1)
    for y1MED = 1:arrSize(2)
        pointCord = [feature1Vals(x1MED); feature2Vals(y1MED)];
        distanceA = sum((muA-pointCord).^2).^0.5;
        distanceB = sum((muB-pointCord).^2).^0.5;
        
        if abs(distanceA - distanceB) < .01
            boundary1 = [boundary1, pointCord];
        end
    end
end 
figure(clusters1)
plot(boundary1(1,:), boundary1(2,:))
polyX = [minValue(1) boundary1(1,length(boundary1)) boundary1(1,1) minValue(1)];
polyY = [minValue(2) boundary1(2,length(boundary1)) boundary1(2,1) minValue(2)];

aInPoly = inpolygon(classA(:,1), classA(:,2), polyX, polyY);
aInBound = numel(classA(aInPoly));

%% MED For Clusters 2
compositeVec = [classC; classD; classE];
maxValue = ceil(max(compositeVec));
minValue = floor(min(compositeVec));

gx = [(muD - muE)', .5*(muE'*muE - muD'*muD)];
vec = [-gx(2)/gx(1) -gx(3)/gx(1)];
refline(vec(1), vec(2))

feature1Vals = minValue(1):0.01:maxValue(1);
feature2Vals = minValue(2):0.01:maxValue(2);

arrSize = [size(feature1Vals,2) size(feature2Vals,2)];
boundary2EC = [];
boundary2DC = [];
boundary2ED = [];

for x2MED = 1:arrSize(1)
    for y2MED = 1:arrSize(2)
        pointCord = [feature1Vals(x2MED); feature2Vals(y2MED)];
        distanceC = sum((muC-pointCord).^2).^0.5;
        distanceD = sum((muD-pointCord).^2).^0.5;
        distanceE = sum((muE-pointCord).^2).^0.5;
        
        if (abs(distanceE - distanceC) < .001) && ((distanceD > distanceE) && (distanceD > distanceC))
            boundary2EC = [boundary2EC, pointCord]; %#ok<AGROW>
        end
        if (abs(distanceD - distanceC) < .001) && ((distanceE > distanceD) && (distanceE > distanceC)) %#ok<*BDSCI>
            boundary2DC = [boundary2DC, pointCord]; %#ok<AGROW>
        end
        if (abs(distanceE - distanceD) < .001) && ((distanceC > distanceE) && (distanceC > distanceD)) 
           boundary2ED = [boundary2ED, pointCord]; %#ok<AGROW>
        end
    end
end 
figure(clusters2)
plot(boundary2ED(1,:),boundary2ED(2,:))
plot(boundary2DC(1,:),boundary2DC(2,:))
plot(boundary2EC(1,:),boundary2EC(2,:))


%% Nearest Neighbour
[X_nn1, Y_nn1, classifier_nn1] = nearestNeighbourFilter(1,classA, 'Class A', classB, 'Class B');
[X_nn2, Y_nn2, classifier_nn2] = nearestNeighbourFilter(1,classC,'Class C',classD,'Class D', classE,'Class E');

% Create test data
classA_test = repmat(muA',[nA, 1])  + randn(nA,2)*chol(sigmaA);
classB_test = repmat(muB',[nB, 1]) + randn(nB,2)*chol(sigmaB);
classC_test = repmat(muC',[nC, 1])  + randn(nC,2)*chol(sigmaC);
classD_test = repmat(muD',[nD, 1]) + randn(nD,2)*chol(sigmaD);
classE_test = repmat(muE',[nE, 1])  + randn(nE,2)*chol(sigmaE);

% Calculate error
nn1_classify = classifyPoints(X_nn1, Y_nn1, classifier_nn1, classA_test, 1, classB_test, 2);
conf_nn1 = confusionmat(nn1_classify(:,1),nn1_classify(:,2));
error_nn1 = size(find(nn1_classify(:,1) ~= nn1_classify(:,2)),1)/size(nn1_classify,1);

nn2_classify = classifyPoints(X_nn2, Y_nn2, classifier_nn2, classC_test, 1, classD_test, 2, classE_test, 3);
conf_nn2 = confusionmat(nn2_classify(:,1),nn2_classify(:,2));
error_nn2 = size(find(nn2_classify(:,1) ~= nn2_classify(:,2)),1)/size(nn2_classify,1);


%% 5 Nearest Neighbour
[X_nn5_1, Y_nn5_1, classifier_nn5_1] = nearestNeighbourFilter(5,classA, 'Class A', classB, 'Class B');
[X_nn5_2, Y_nn5_2, classifier_nn5_2] = nearestNeighbourFilter(5,classC,'Class C',classD,'Class D', classE,'Class E');

% Calculate error
nn5_1_classify = classifyPoints(X_nn5_1, Y_nn5_1, classifier_nn5_1, classA_test, 1, classB_test, 2);
conf_nn5_1 = confusionmat(nn5_1_classify(:,1),nn5_1_classify(:,2));
error_nn5_1 = size(find(nn5_1_classify(:,1) ~= nn5_1_classify(:,2)),1)/size(nn5_1_classify,1);

nn5_2_classify = classifyPoints(X_nn5_2, Y_nn5_2, classifier_nn5_2, classC_test, 1, classD_test, 2, classE_test, 3);
conf_nn5_2 = confusionmat(nn5_2_classify(:,1),nn5_2_classify(:,2));
error_nn5_2 = size(find(nn5_2_classify(:,1) ~= nn5_2_classify(:,2)),1)/size(nn5_2_classify,1);
