%% Aliasing Demonstration
% Senan Ebrahim
% August 18, 2016
% HST 015

% Run the script to see an example of spatial aliasing. This is why we
% prefer to smooth the data than to simply downsample (see anti-aliasing filters).

% Creating pic 1 = original
stan1 = imread('Robbins1.png');

% Creating pic 2 = different signal than pic 1 (adding black bars)
indices = 2:2:size(stan1,2);
z = zeros(size(stan1,1),size(stan1,2));
stanbars=stan1;
stanbars(:,indices) = z(:,indices);

% Here we create downsampled versions of both pics (taking every other
% datapoint, you can change 2 to modify downsampling factor) 
standown = stan1(1:2:end,1:2:end);
stanbarsdown = stanbars(1:2:end,1:2:end);

% Here are the original pics, note that they appear different 
% due to the bars  
figure;
imshow([stan1,stanbars]);
title('original pic 1                |                original pic 2');

% Here are the downsampled pics, now they appear the same
figure;
imshow([standown,stanbarsdown])
title('downsampled pic 1   |   downsampled pic 2');

% The loss of information such that we cannot distinguish the two pictures
% is aliasing. We generally do not want to have this happen (assuming the
% difference between the pictures is important to us).