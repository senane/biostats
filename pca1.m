
function [pc, scores, vars, exps] = pca1(dat)


% Senan Ebrahim
% August 18, 2016
% HST 015

% Based on Jonathon Shlens PCA Tutorial Appendix B

% pc is the principal component matrix
% scores is our projection of the original data
% vars is the variances of the PCs (equivalent to eigenvalues)
% exps is the percent variance explained by each PC


% data input is in the form (n x p), with n trials of p dimensions
%data = data';

[n, p] = size(dat);
dat = dat - repmat(mean(dat),n,1);

covar = cov(dat);

[pc,vars] = eig(covar);

%need to extract the variances of the PCs (rest of matrix is 0 since PCs
%are orthonormal)
vars = diag(vars);

[noop, rindices] = sort(-1*vars);
vars = vars(rindices);
pc = pc(:,rindices);

exps = 100*vars/sum(vars);

scores = dat*pc;
end