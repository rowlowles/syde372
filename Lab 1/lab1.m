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
bInPoly = inpolygon(classB(:,1), classB(:,2), polyX, polyY);
bInBound = nB - numel(classB(bInPoly));


%% MED For Clusters 2
compositeVec = [classC; classD; classE];
maxValue = ceil(max(compositeVec));
minValue = floor(min(compositeVec));

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