function [feature1Vals, feature2Vals, MLClassificationCDE] = ML(classA, meanC, sigmaC, classD, meanD, sigmaD, classE, meanE, sigmaE)
% set-up grid in featurespace to be populated with classifications
compositeVec = [classA; classD; classE];
maxValue = ceil(max(compositeVec));
minValue = floor(min(compositeVec));
feature1Vals = minValue(1):1:maxValue(1);
feature2Vals = minValue(2):1:maxValue(2);
arrSize = [size(feature2Vals,2) size(feature1Vals,2)];
MLClassificationCDE = zeros(arrSize);

% Set the a priori class probabilities proportional to the number of 
% samples in each class.
priorC = 1/3;
priorD = 1/3;
priorE = 1/3;

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
            MLClassificationCDE(j,i) = 1; %class C
        elseif(logPosteriorD(j,i) > logPosteriorC(j,i) && logPosteriorD(j,i) > logPosteriorE(j,i))
            MLClassificationCDE(j,i) = 2; %class D
        elseif(logPosteriorE(j,i) > logPosteriorC(j,i) && logPosteriorE(j,i) > logPosteriorD(j,i))
            MLClassificationCDE(j,i) = 3; %class E
        else
            MLClassificationCDE(j,i) = 0; %boundary!
        end
    end
end

end

