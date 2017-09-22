%% PCA: An Overview
% Senan Ebrahim
% August 17, 2017
% HST 015

% Goal 1 today: take a deeper look under the hood of matrix manipulation in
% matlab.

% Goal 2: get a rough idea of PCA as a versatile data analysis tool.

% First, we will cover one implementation of PCA. Then we will move on to apply PCA to a few data sets.

%% Covariance
% Covariance is defined as Cov(X,Y) = E(X-E(X))*E(Y-E(Y))
x = rand(3,1)
y = rand(3,1)

covariance = sum((x-mean(x)).*(y-mean(y)))/(length(x)-1)

covar = @(u,v) sum(((u-mean(u)).*(v-mean(v))))/(length(u)-1);

mcov = cov(x,y);

scov = covar(x,y);

covariance == scov && scov == mcov(1,2)

%% Eigenvectors

A = [1 2; 3 4];
l = eig(A)
[v, d] = eig(A)

A2 = int8(v*d*inv(v)) % If v is eigenvectors, we should be able to reconstitute our dataset

A == A2 
% It is now true but if you use an old version of MATLAB comes up
% false because of floating point error.

isequal(A,A2)

%% PCA implementation using eigenvectors
% See function pca1



%% PCA on Multivariate CHD Dataset to test implementation
% Dataset modified from UCI ML Database
load('PCAdata.mat')

c1 = pca1(chd);
c2 = pca(chd,'Algorithm','eig');
c3 = pca(chd,'Algorithm','svd');

test = isequal(c1,c2)

% Test what the biggest difference is (order of magnitude for norm)
delts = abs(c1)-abs(c2);
max(delts(:))
norm(delts)

% Returns the Euclidian norm (p=2) which is the magnitude of each vector in
% the delta space
% Creating an inline function so we do not have to keep retyping
normabsd = @(x,y) norm(abs(x)-abs(y));
normabsd(c1,c2)
normabsd(c2,c3)
normabsd(c1,c3)

% Is our magnitude less than 1e-12? If so we have effectively the same answer since we
% are at the limit of double floating point precision.


%% PCA on Bivariate Retina Dataset
% Dataset modified from UCI ML Database
% Data is measurements of optic disc diameter and optic disc to macula distance 200 patients with diabetic retinopathy

% Visualize data before transform

% Ensure normality
figure;
subplot(2,1,1);
histogram(retina(:,1));
xlabel(retina_labels{1});
ylabel('freq');
subplot(2,1,2);
histogram(retina(:,2));
xlabel(retina_labels{2});
ylabel('freq');

%View scatter
figure;
subplot(2,2,1);
scatter(retina(:,1),retina(:,2), 'x');
xlabel(retina_labels{1});
ylabel(retina_labels{2});

% Run PCA
[Rpc, Rscores, Rvars, Rexps] = pca1(retina);

% Visualize data after transform
subplot(2,2,2);
scatter(Rscores(:,1),Rscores(:,2),'x');
set(gca,'XLim',[-0.3 0.3],'YLim',[-0.1 0.1]);
xlabel('PC1');
ylabel('PC2');

% Visualize data collapsed onto a single PC
subplot(2,2,3);
scatter(Rscores(:,1),zeros(size(Rscores,1),1),'x');
set(gca,'XLim',[-0.3 0.3],'YLim',[-0.1 0.1]);
xlabel('PC1');

% Visualize data rotated back onto original axes
Rv1 = Rpc(:,1);
retinared = retina*Rv1; %reduced data set 
retinadec = retinared*Rv1'; %decompressed data set
subplot(2,2,4);
scatter(retinadec(:,1),retinadec(:,2), 'x');
xlabel(retina_labels{1});
ylabel(retina_labels{2});

% Visualize covariance of variables and covariance of PCs
figure;
imagesc(cov(retina));
colorbar;
title('Covariance of Original Retina Data')

set(gca,'XLim',[0.5 2.5],'XTick',0.5:0.5:2.5,'YLim',[0.5 2.5],'YTick',0.5:0.5:2.5,...
    'XTickLabel',{'',retina_labels{1},'',retina_labels{2},''},...
    'YTickLabel',{'',retina_labels{1},'',retina_labels{2},''});

figure;
imagesc(cov(Rscores));
colorbar;
title('Covariance of Retina PCs')
xlabel('PC 1 ID #')
ylabel('PC 2 ID #')


Rexps % this shows us the percent of variance explained by the 2 PCs

% Q: Is PCA appropriate on this dataset?
% A: On the one hand, our eigenvalues are way higher for PC 1 for PC 2 such
% that PC 1 accounts for 99% of the variability.
% On the other, the variation in the ratio is actually pretty important for 
% diagnosing disease severity in diabetic retinopathy.
% 


%% PCA on CHD Dataset Analysis

% In-class exercise - go to announcements
%
% 
% Visualize covariance and correlation matrices of variables 
figure;
imagesc(cov(chd));
colorbar;
title('Covariance of CHD Data')


figure;
imagesc(corr(chd));
colorbar;
title('Correlation of CHD Data')
%% 
% Run PCA on both the original data and the original data normalized

% PCA on original data set
[Cpc, Cscores, Cvars, Cexps] = pca1(chd);

figure
plot(Cexps, '-o')
title('Scree Plot of PC Contributions (Cov)');
xlabel('PC ID');
ylabel('Percent Variance Explained');

% Percent variance explained by the first three PCs
first3expcov=sum(Cexps(1:3))

% PCA on on normalized data set
[Cpc1, Cscores1, Cvars1, Cexps1] = pca1(zscore(chd));

figure;
plot(Cexps1, '-o');
title('Scree Plot of PC Contributions (Corr)');
xlabel('PC ID');
ylabel('Percent Variance Explained');

% Percent variance explained by the first three PCs
first3expcorr=sum(Cexps1(1:3))

% Q: Was PCA a useful technique for this dataset? Which one is 'better'?
% Q: What does each PC represent?