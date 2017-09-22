%% Machine Learning Lecture
% Senan Ebrahim
% August 31, 2017
% HST 015


%% Regression of cancer data

% Our goal is to predict survival time based on tumor size and age.


load('MLdata.mat');
rng('default');

[betas deviance stats] = glmfit(cancer(:,1:2),cancer(:,4));
predsurvival = glmval(betas,cancer(:,1:2),'identity');

plot3(cancer(:,1),cancer(:,2),predsurvival, 'ro')
hold on
plot3(cancer(:,1),cancer(:,2),cancer(:,4), 'gx')
shg
xlabel('Age')
ylabel('Tumor Size')
zlabel('Survival Time')

% Can see a plane there, could do pcfitplane if you want it

% BACK TO SLIDES

% We can assess model performance by looking at deviance or sum of squares 
RSS = sum(stats.resid.^2);
deviance == RSS

% Let us validate our regression by looking at the residuals
figure;
scatter3(cancer(:,1),cancer(:,2),stats.resid);
xlabel('Age')
ylabel('Tumor Size')
zlabel('Residuals')

figure;
histogram(stats.resid)
title('Histogram of Residuals')

% Various statistical methods to do this rigorously

% Best way to validate your model is really cross-validation, we will use
% k-fold
predmodel = fitrlinear(cancer(:,1:2),cancer(:,4),'Learner','leastsquares','CrossVal','on','KFold',10)
predmodel1 = predmodel.Trained{1};

predsurvival = kfoldPredict(predmodel);

%MSE calculation
predsurvivalloss = kfoldLoss(predmodel)

% Matlab function with built-in regularization: lassoglm
% But here we will continue using fitrlinear because it has built-in
% cross-validation for our lasso

% BACK TO SLIDES

% Test a single lambda
Lambda = 0.1;
predmodelreg = fitrlinear(cancer(:,1:2),cancer(:,4),'KFold',5,'Lambda',Lambda, 'Learner','leastsquares','Regularization','lasso');
predsurvivalreg = kfoldPredict(predmodelreg);
predsurvivalregloss = kfoldLoss(predmodelreg)

% Initialize some lambdas evenly spaced log-scale to test
Lambdas = logspace(-5,5,20);
predmodelregs = fitrlinear(cancer(:,1:2),cancer(:,4),'KFold',5,'Lambda',Lambdas, 'Learner','leastsquares','Regularization','lasso');
predsurvivalregslosses = kfoldLoss(predmodelregs)

% Identify best-performing lambda
scatter(log(Lambdas),predsurvivalregslosses);
xlabel('Log(Lambda)');
ylabel('Loss (MSE)');
[ val ind ] = min(predsurvivalregslosses)
Lambdas(ind)


% use App

% BACK TO SLIDES



%% Classification of cancer data

% Our goal is to predict malignancy status based on tumor size and age.
rng('default');

% We can create a single decision tree that does this
maligclassifier = fitctree(cancer(:,1:2),cancer(:,3));

% Visualize the tree
view(maligclassifier,'Mode','graph')

% Evaluate our model using ROC / AUC
[~, score] = resubPredict(maligclassifier);
[ROCx, ROCy, ROCt, AUC] = perfcurve(cancer(:,3),score(:,2),1);
figure
plot(ROCx,ROCy);
AUC
title('AUC')

% Now we try predicting malignancy status using a random forest of ten
% trees
maligclassifierbagged = TreeBagger(10,cancer(:,1:2),cancer(:,3),'OOBPrediction','On','Method','classification');

% Visualize one tree - the first one
view(maligclassifierbagged.Trees{1},'Mode','graph')

% Prediction is aggregated across all trees so let's see how performance
% improves with number of trees
figure;
oobErrorBaggedEnsemble = oobError(maligclassifierbagged);
plot(oobErrorBaggedEnsemble)
xlabel('Number of grown trees');
ylabel('Out-of-bag classification error');


% use App


%% In Class Exercise
% 
% Predict survival time given age, tumor size, and malignancy status
rng('default');

[betasm deviancem statsm] = glmfit(cancer(:,1:3),cancer(:,4));
predsurvivalm = glmval(betasm,cancer(:,1:3),'identity');

% You could do a PCA on this and find fewer axes but for now let's just do
% pairwise scatters to visualize residuals

figure;
subplot(221)
scatter(cancer(:,1),statsm.resid);
xlabel('Age')
ylabel('Residuals')

subplot(222)
scatter(cancer(:,2),statsm.resid);
xlabel('Tumor Size')
ylabel('Residuals')

subplot(223)
scatter(cancer(:,3),statsm.resid);
xlabel('Malignancy Status')
ylabel('Residuals')

subplot(224)
histogram(statsm.resid)
title('Histogram of Residuals')

% 5-fold cross-validate the model
% Best way to validate your model is called cross-validation (to slides)
predmodelm = fitrlinear(cancer(:,1:2),cancer(:,4),'Learner','leastsquares','CrossVal','on','KFold',5)
predmodelm1 = predmodelm.Trained{1};

predsurvival = kfoldPredict(predmodelm);

%MSE calculation
predsurvivalloss = kfoldLoss(predmodelm)

% Is it overfit? We do not know until we test on outside data.

% Use a lasso or a ridge regression to regularize
% Lasso
Lambda = 0.1;
predmodellasso = fitrlinear(cancer(:,1:2),cancer(:,4),'KFold',5,'Lambda',Lambda, 'Learner','leastsquares','Regularization','lasso');
predsurvivallasso = kfoldPredict(predmodellasso);
predsurvivallasso = kfoldLoss(predmodellasso)

% Ridge
% Test a single lambda
predmodelridge = fitrlinear(cancer(:,1:2),cancer(:,4),'KFold',5, 'Learner','leastsquares','Regularization','ridge');
predsurvivalridge = kfoldPredict(predmodelridge);
predsurvivalridge = kfoldLoss(predmodelridge)

% Reframe survival time as a classification problem of whether or not the patient survives 10+ years (binarize)
survivaltime = cancer(:,4)>10*365;

% Build a random forest classifier and assess its performance.
survivalclassifierbagged = TreeBagger(50,cancer(:,1:3),survivaltime,'OOBPrediction','On','Method','classification');

% Visualize one tree - the first one
view(survivalclassifierbagged.Trees{1},'Mode','graph')

% Prediction is aggregated across all trees so let's see how performance
% improves with number of trees
figure;
surverror = oobError(survivalclassifierbagged);
plot(surverror)
xlabel('Number of grown trees');
ylabel('Out-of-bag classification error');


%% Clustering of genomic data
% Here we have gene expression levels for 614 genes from 7 individuals

rng('default')

figure
[~,score,~,~,explainedVar] = pca(genomic);
plot(explainedVar, '-o')
title('PC Contributions to Variance')
ylabel('Percent Variance Explained')
xlabel('PC ID')

% Retain first two principal components
PCscore = score(:,1:2);

% Visualize data in two PCs
figure
scatter(PCscore(:,1),PCscore(:,2));
xlabel('PC 1');
ylabel('PC 2');
title('Genetic Expression Data in PC Space');

% Perform clustering
nclusters = 2;
[clusters, centroid] = kmeans(PCscore,nclusters);

% Do with 2, then with 10

% Are they equal sizes? Not in number of data points! Looking at euclidian
% distance.

% Now plot with cluster grouping
gscatter(PCscore(:,1),PCscore(:,2),clusters);
legend('location','southeast')
xlabel('First Principal Component');
ylabel('Second Principal Component');
title('Genetic Expression Data in PC Space with Clustering');

%% Clustering a histological image

nissl = imread('Nissl.JPG');
figure
imshow(nissl)
title('original image');

% Our eye sees two colors here


% Now we need to think about what feature space to use to describe the
% color and brightness independently

% By default our image is described in RGB format

% There is a color space called CIELAB = L*a*b that separates out
% luminosity/brightness, red-green value, and blue-yellow value

cform = makecform('srgb2lab');
labnissl = applycform(nissl,cform);
ab = double(labnissl(:,:,2:3));
nrows = size(ab,1);
ncols = size(ab,2);
ab = reshape(ab,nrows*ncols,2);

% Since we are looking for two colors
nColors = 2;

% Perform k-means clustering
[cluster_idx cluster_center] = kmeans(ab,nColors);

pixel_labels = reshape(cluster_idx,nrows,ncols);
figure
imshow(pixel_labels,[]);
title('image labeled by cluster index');
shg
segmented_images = cell(1,3);
rgb_label = repmat(pixel_labels,[1 1 3]);
for k = 1:nColors
color = nissl;
color(rgb_label ~= k) = 0;
segmented_images{k} = color;
end
figure
imshow(segmented_images{1});
title('cluster 1');
figure
imshow(segmented_images{2});
title('cluster 2');


