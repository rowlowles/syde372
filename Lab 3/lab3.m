close all;

load feat.mat
featureSpace = figure;

subplot(1,3,1);
title("Feature Space for 2x2 image blocks");
xlabel("Feature 1");
ylabel("Feature 2");
aplot(f2);

subplot(1,3,2);
title("Feature Space for 8x8 image blocks");
xlabel("Feature 1");
ylabel("Feature 2");
aplot(f8);

subplot(1,3,3);
title("Feature Space for 32x32 image blocks");
xlabel("Feature 1");
ylabel("Feature 2");
aplot(f32);

%% Create classifiers
imageClassification(f2,f2t,"f2","Image Classification Using 2x2 Image Blocks");
imageClassification(f8,f8t,"f8","Image Classification Using 8x8 Image Blocks");
imageClassification(f32,f32t,"f32","Image Classification Using 32x32 Image Blocks");
