close all;

load feat.mat
% featureSpace = figure;

% subplot(1,3,1);
figure
title("Feature Space for 2x2 image blocks");
xlabel("Feature 1");
ylabel("Feature 2");
aplot(f2);

% subplot(1,3,2);
figure
title("Feature Space for 8x8 image blocks");
xlabel("Feature 1");
ylabel("Feature 2");
aplot(f8);

% subplot(1,3,3);
figure
title("Feature Space for 32x32 image blocks");
xlabel("Feature 1");
ylabel("Feature 2");
aplot(f32);

%% Create classifiers
[classifier_f2, X_f2, Y_f2] = imageClassification(f2,f2t,"f2","Image Classification Using 2x2 Image Blocks");
[classifier_f8, X_f8, Y_f8] = imageClassification(f8,f8t,"f8","Image Classification Using 8x8 Image Blocks");
[classifier_f32, X_f32, Y_f32] = imageClassification(f32,f32t,"f32","Image Classification Using 32x32 Image Blocks");

%%  Image Classification and Segmentation
cimage = zeros(256,256);

for i = 1:256
    for j = 1:256
        cimage(i,j) = findClassification(X_f8,Y_f8,classifier_f8,...
            multf8(i,j,1), multf8(i,j,2));
    end
end

im_cl_sg_fig = figure;
subplot(1,2,1);
imagesc(cimage);
title("Classification of Composite Image of Multiple Textures");
xlabel("Horizontal Pixel");
ylabel("Vertical Pixel");
colorbar();

subplot(1,2,2);
imshow(multim, []);
title("Composite Image of Multiple Textures");
xlabel("Horizontal Pixel");
ylabel("Vertical Pixel");
axis on;

saveas(im_cl_sg_fig, "comp_img.png");

% Abstract class image
numOfTicks = 9;
xtick_num = linspace(0.0035, 0.2042, numOfTicks);
ytick_num = linspace(0.0062, 0.2036, numOfTicks);

xtick = cell(numOfTicks-2, 1);
ytick = cell(numOfTicks-2, 1);

for i = 2:(numOfTicks-1) 
    xtick{i-1} = xtick_num(i);
    ytick{i-1} = ytick_num(i);
end

abstract_art = figure;
imagesc(classifier_f8);
ax = gca;
ax.YDir = 'normal';
title("Classifier for 8x8 Blocks of Texture Images");
xlabel("Feature 1");
ylabel("Feature 2");
xticklabels(xtick);
yticklabels(ytick);
colorbar();

saveas(abstract_art, "abstract_art.png");

%% K Means
kmeans(10,f32(1:2,:)');