%% PCA on PBC Data
%
% Data Courtesy of Mayo Clinic (modified)
% HST 190
% SE 8/17/2016


%% Close figures and clear workspace
% !!! Make sure you have any figures and workspace variables you need saved
% before running the script!
clear;      
close all;  

%% Load data
datatable=readtable('PBC.xlsx');
data=table2array(datatable);
varlist=datatable.Properties.VariableNames;


%% Calculate and plot covariance matrix (Figure 1)
% For part (b) change this line so that we visualize the correlations
% rather than the covariances
dataCov = cov(data); %YOUR CODE HERE

figure;
imagesc(dataCov);
colorbar;

% Change Title to match plot
title('Covariance of 12 Variables in Primary Biliary Cirrhosis'); %YOUR CODE HERE
xlabel('Variable 1 ID');
ylabel('Variable 2 ID');
set(gca,'XLim',[0 12],'XTick',1:12,'YLim',[0 12],'YTick',1:12,'XTickLabel',varlist, 'YTickLabel',varlist);

% Making figure nice
set(gcf, 'Position',[5 5 960 960 ],'PaperUnits','inches','PaperSize',[8 8],'Color', [1 1 1]);                  
set(findall(gcf,'-property','FontSize'),'FontSize',12);
set(findall(gcf,'-property','Box'),'Box','off');

%% Run PCA
% For part (c) change this line so that the PCA is done using the
% correlations instead of covariances (there are multiple ways to do this)
[coeff,score,latent,tsquared,lambdas,mu] = pca(data); %YOUR CODE HERE

%% Plot Eigenvalues (Figure 2)
figure;
plot(lambdas, '-o');
title('Scree Plot of PC Contributions');
xlabel('PC ID');
ylabel('Eigenvalue');

%% Find and plot the covariance matrix of the PCs (Figure 3)
% Fill in the code for pcCov by calculating the covariance of the PCs
pcCov = []; %YOUR CODE HERE
figure;
imagesc(pcCov);
colorbar;
title('Covariance of 12 PCs');
xlabel('PC 1 ID');
ylabel('PC 2 ID');

% Making figure nice
set(gcf, 'Position',[5 5 960 960 ],'PaperUnits','inches','PaperSize',[8 8],'Color', [1 1 1]);                  
set(findall(gcf,'-property','FontSize'),'FontSize',12);
set(findall(gcf,'-property','Box'),'Box','off');
