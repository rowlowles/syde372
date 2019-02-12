% Lab 1
close all;
clear all;
% Generate clusters
% Class A
nA = 200;
muA = [5;10];
sigmaA = [8 0; 0 4];
classA = repmat(muA',[nA, 1]) + randn(nA,2)*chol(sigmaA);
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

%% MED For Clusters 1
minValue = floor(min(min(classA, classB)));
maxValue = ceil(max(max(classA, classB)));

feature1Vals = minValue(1):0.1:maxValue(1);
feature2Vals = minValue(2):0.1:maxValue(2);

[X_MED_1, Y_MED_1] = meshgrid(feature1Vals, feature2Vals);

arrSize = [size(feature2Vals,2) size(feature1Vals,2)];

classifier_MED_1 = zeros(arrSize);
boundary1 = [];

for x1MED = 1:arrSize(2)
    for y1MED = 1:arrSize(1)
        pointCord = [feature1Vals(x1MED); feature2Vals(y1MED)];
        distanceA = sum((muA-pointCord).^2).^0.5;
        distanceB = sum((muB-pointCord).^2).^0.5;
        
        if abs(distanceA - distanceB) < .01
            boundary1 = [boundary1, pointCord];
        end
        
        if distanceA < distanceB
            classifier_MED_1(y1MED,x1MED) = 1;
        else
            classifier_MED_1(y1MED,x1MED) = 2;
        end
    end
end

figure;
plotClasses(classA,'Class A',classB,'Class B');
hold on;
plot(boundary1(1,:), boundary1(2,:))
polyX = [minValue(1) boundary1(1,length(boundary1)) boundary1(1,1) minValue(1)];
polyY = [minValue(2) boundary1(2,length(boundary1)) boundary1(2,1) minValue(2)];

oneMat = ones(200,1);

aInPoly = double(inpolygon(classA(:,1), classA(:,2), polyX, polyY));
bInPoly = (oneMat - inpolygon(classB(:,1), classB(:,2), polyX, polyY));

aInPoly(aInPoly == 0) = [2]
bInPoly = bInPoly+1;

attachedMat = [oneMat aInPoly; oneMat+1 bInPoly];
confMatrixA = confusionmat(attachedMat(:,1), attachedMat(:,2))
errorClass = size(find(attachedMat(:,1) ~= attachedMat(:,2)),1)/size(attachedMat,1)

%% MED For Clusters 2
compositeVec = [classC; classD; classE];
maxValue = ceil(max(compositeVec));
minValue = floor(min(compositeVec));

% gx = [(muD - muE)', .5*(muE'*muE - muD'*muD)];
% vec = [-gx(2)/gx(1) -gx(3)/gx(1)];
% refline(vec(1), vec(2))

feature1Vals = minValue(1):0.01:maxValue(1);
feature2Vals = minValue(2):0.01:maxValue(2);

arrSize = [size(feature1Vals,2) size(feature2Vals,2)];
classifierMat = zeros(arrSize);

boundary2EC = [];
boundary2DC = [];
boundary2ED = [];


[X_MED, Y_MED] = meshgrid(feature1Vals, feature2Vals);
classifier_MED_2 = zeros(arrSize);

for x2MED = 1:arrSize(2)
    for y2MED = 1:arrSize(1)
        pointCord = [feature1Vals(x2MED); feature2Vals(y2MED)];
        distanceC = sum((muC-pointCord).^2).^0.5;
        distanceD = sum((muD-pointCord).^2).^0.5;
        distanceE = sum((muE-pointCord).^2).^0.5;
        
        if (distanceC < distanceE) && (distanceC < distanceD)
            classifierMat(x2MED, y2MED) = 1;
        end
        if (distanceD < distanceE) && (distanceD < distanceC)
            classifierMat(x2MED, y2MED) = 2;
        end
        if (distanceE < distanceC) && (distanceE < distanceD)
            classifierMat(x2MED, y2MED) = 3;
        end
        
        if (abs(distanceE - distanceC) < .001) && ((distanceD > distanceE) && (distanceD > distanceC))
            boundary2EC = [boundary2EC, pointCord]; %#ok<AGROW>
        end
        if (abs(distanceD - distanceC) < .001) && ((distanceE > distanceD) && (distanceE > distanceC)) %#ok<*BDSCI>
            boundary2DC = [boundary2DC, pointCord]; %#ok<AGROW>
        end
        if (abs(distanceE - distanceD) < .001) && ((distanceC > distanceE) && (distanceC > distanceD)) 
           boundary2ED = [boundary2ED, pointCord]; %#ok<AGROW>
        end
        
        if ((distanceC < distanceD) && (distanceC < distanceE))
            classifier_MED_2(y2MED,x2MED) = 1;
        elseif((distanceD < distanceE)) 
             classifier_MED_2(y2MED,x2MED) = 2;   
        else
             classifier_MED_2(y2MED,x2MED) = 3;
        end
        
    end
end

classifiedPoints = classifyPoints(X_MED, Y_MED, classifierMat, classC, 1, classD, 2, classE, 3);
conf_MED = confusionmat(classifiedPoints(:,1),classifiedPoints(:,2))
error_MED = size(find(classifiedPoints(:,1) ~= classifiedPoints(:,2)),1)/size(classifiedPoints,1)

figure;
plotClasses(classC,'Class C',classD,'Class D', classE,'Class E');
hold on; 
plot(boundary2ED(1,:),boundary2ED(2,:)) % Bottom Right
plot(boundary2DC(1,:),boundary2DC(2,:)) % Top line
plot(boundary2EC(1,:),boundary2EC(2,:)) % Bottom left

%% GED Classifier
[X_ged1, Y_ged1, classifier_ged1] = GEDFilter(classA, muA, sigmaA, 'Class A', classB, muB, sigmaB, 'Class B');
[X_ged2, Y_ged2, classifier_ged2] = GEDFilter(classC, muC, sigmaC, 'Class C', classD, muD, sigmaD, 'Class D', classE, muE, sigmaE, 'Class E');

% Calculate error
ged1_classify = classifyPoints(X_ged1, Y_ged1, classifier_ged1, classA, 1, classB, 2);
conf_ged1 = confusionmat(ged1_classify(:,1), ged1_classify(:,2));
error_ged1 = size(find(ged1_classify(:,1) ~= ged1_classify(:,2)),1)/size(ged1_classify,1);

ged2_classify = classifyPoints(X_ged2, Y_ged2, classifier_ged2, classC, 1, classD, 2, classE, 3);
conf_ged2 = confusionmat(ged2_classify(:,1), ged2_classify(:,2));
error_ged2 = size(find(ged2_classify(:,1) ~= ged2_classify(:,2)),1)/size(ged2_classify,1);

%% MAP Classifier
% Maximum A Posterioi (MAP), using the true statistics. Set the a
% priori class probabilities proportional to the number of samples in
% each class.

% close all

% Classes A & B

% set-up grid in featurespace to be populated with classifications
minValue = floor(min(min(classA, classB)));
maxValue = ceil(max(max(classA, classB)));
feature1Vals = minValue(1):0.1:maxValue(1);
feature2Vals = minValue(2):0.1:maxValue(2);
arrSize = [size(feature2Vals,2) size(feature1Vals,2)];
MAPClassificationAB = zeros(arrSize);


% Set the a priori class probabilities proportional to the number of 
% samples in each class.
priorA = nA / (nA + nB);
priorB = nB / (nA + nB);

%initialize posteriors
logPosteriorA = zeros(arrSize);
logPosteriorB = zeros(arrSize);

for i = 1:size(feature1Vals,2)
    for j = 1:size(feature2Vals,2)
        x = [feature1Vals(i), feature2Vals(j)];
        logLikihoodA = -log(2*pi*sqrt(det(sigmaA))) -0.5 * (x - meanA) * inv(sigmaA) * (x - meanA)';
        logLikihoodB = -log(2*pi*sqrt(det(sigmaB))) -0.5 * (x - meanB) * inv(sigmaB) * (x - meanB)';
        logPosteriorA(j,i) = logLikihoodA + log(priorA);
        logPosteriorB(j,i) = logLikihoodB + log(priorB);
        
        if(logPosteriorA(j,i) > logPosteriorB(j,i))
            MAPClassificationAB(j,i) = 1; %class A            
        elseif(logPosteriorA(j,i) < logPosteriorB(j,i))
            MAPClassificationAB(j,i) = 2; %class B
        else
            MAPClassificationAB(j,i) = 0; %boundary!
        end
        
    end
end

%plotting
figure
hold on;
set(gca, 'ydir', 'normal');
%axis equal;
plotClasses(classA,'Class A',classB,'Class B');
plotStdContours([1], meanA, sigmaA, meanB, sigmaB);
[X_MAP_1, Y_MAP_1] = meshgrid(feature1Vals, feature2Vals);
contour(X_MAP_1,Y_MAP_1,MAPClassificationAB,'DisplayName','MAP boundary')

% Classes C, D, E

% set-up grid in featurespace to be populated with classifications
compositeVec = [classA; classB; classC; classD; classE];
maxValue = ceil(max(compositeVec));
minValue = floor(min(compositeVec));
feature1Vals = minValue(1):0.1:maxValue(1);
feature2Vals = minValue(2):0.1:maxValue(2);
arrSize = [size(feature2Vals,2) size(feature1Vals,2)];
MAPClassificationCDE = zeros(arrSize);

% Set the a priori class probabilities proportional to the number of 
% samples in each class.
priorC = nC / (nC + nD + nE);
priorD = nD / (nC + nD + nE);
priorE = nE / (nC + nD + nE);

%initialize posteriors
logPosteriorC = zeros(arrSize);
logPosteriorD = zeros(arrSize);
logPosteriorE = zeros(arrSize);

for i = 1:size(feature1Vals,2)
    for j = 1:size(feature2Vals,2)
        x = [feature1Vals(i), feature2Vals(j)];
        
        logLikihoodC = -log(2*pi*sqrt(det(sigmaC))) -0.5 * (x - meanC) * inv(sigmaC) * (x - meanC)';
        logLikihoodD = -log(2*pi*sqrt(det(sigmaD))) -0.5 * (x - meanD) * inv(sigmaD) * (x - meanD)';
        logLikihoodE = -log(2*pi*sqrt(det(sigmaE))) -0.5 * (x - meanE) * inv(sigmaE) * (x - meanE)';
        
        logPosteriorC(j,i) = logLikihoodC + log(priorC);
        logPosteriorD(j,i) = logLikihoodD + log(priorD);
        logPosteriorE(j,i) = logLikihoodE + log(priorE);
        
        if(logPosteriorC(j,i) > logPosteriorD(j,i) && logPosteriorC(j,i) > logPosteriorE(j,i))
            MAPClassificationCDE(j,i) = 1; %class C
        elseif(logPosteriorD(j,i) > logPosteriorC(j,i) && logPosteriorD(j,i) > logPosteriorE(j,i))
            MAPClassificationCDE(j,i) = 2; %class D
        elseif(logPosteriorE(j,i) > logPosteriorC(j,i) && logPosteriorE(j,i) > logPosteriorD(j,i))
            MAPClassificationCDE(j,i) = 3; %class E
        else
            MAPClassificationCDE(j,i) = 0; %boundary!
        end
    end
end

%plotting
figure
hold on;
set(gca, 'ydir', 'normal');
plotClasses(classC,'Class C',classD,'Class D', classE,'Class E');
plotStdContours([1], meanC, sigmaC, meanD, sigmaD, meanE, sigmaE);
[X_MAP_2, Y_MAP_2] = meshgrid(feature1Vals, feature2Vals);
contour(X_MAP_2,Y_MAP_2,MAPClassificationCDE,'DisplayName','MAP boundary')

% Calculate error
MAP_1_classify = classifyPoints(X_MAP_1, Y_MAP_1, MAPClassificationAB, classA, 1, classB, 2);
conf_MAP_1 = confusionmat(MAP_1_classify(:,1),MAP_1_classify(:,2));
error_MAP_1 = size(find(MAP_1_classify(:,1) ~= MAP_1_classify(:,2)),1)/size(MAP_1_classify,1);

MAP_2_classify = classifyPoints(X_MAP_2, Y_MAP_2, MAPClassificationCDE, classC, 1, classD, 2, classE, 3);
conf_MAP_2 = confusionmat(MAP_2_classify(:,1),MAP_2_classify(:,2));
error_MAP_2 = size(find(MAP_2_classify(:,1) ~= MAP_2_classify(:,2)),1)/size(MAP_2_classify,1);


%% Plot them all together
figure;
plotClasses(classA,'Class A',classB,'Class B');
hold on; 
plotStdContours([1], meanA, sigmaA, meanB, sigmaB);
hold on;
contour(X_MAP_1,Y_MAP_1,MAPClassificationAB,'DisplayName','MAP boundary');
hold on;
contour(X_ged1,Y_ged1,classifier_ged1,'DisplayName','GED boundary');
hold on;
contour(X_MED_1,Y_MED_1,classifier_MED_1,'DisplayName','MED boundary');

figure;
plotClasses(classC,'Class C',classD,'Class D', classE,'Class E');
hold on; 
plotStdContours([1], meanC, sigmaC, meanD, sigmaD, meanE, sigmaE);
hold on;
contour(X_MAP_2,Y_MAP_2,MAPClassificationCDE,'DisplayName','MAP boundary');
hold on;
contour(X_ged2,Y_ged2,classifier_ged2,'DisplayName','GED boundary');
hold on;
contour(X_MED_2,Y_MED_2,classifier_MED_2,'DisplayName','MED boundary');


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
% nn1_classify is a 2d array, where column 1 is a list of points with val
% equalling to which class they belong to in reality. Column 2 is the list
% of points classified according to what the algo thinks they are 
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
