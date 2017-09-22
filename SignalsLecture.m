%% Temporal Signal Analysis
% Senan Ebrahim
% September 1, 2016
% HST 015

% Senan Ebrahim

%% % Time Series Data

% We will start with several examples of time series data (all are signals)
load('Session8.mat');

%% Signal 1: Temperature over time (temps vs. temptimes)
% Temperature (F) every 12 hours for Boston, MA for one year (simulated)
% Time here is in days
plot (temptimes, temps);
title ('Temperature in Boston');
xlabel('Time (days)');
ylabel('Temp (F)');
shg;

%% Signal 2: Protein expression over time
% This is simulated data but for a real experiment done (see Aymoz et al., 2016)
% Expression levels of fluorescent reporter for a gene over 100 minutes

% protexps is an array of four proteins and their expression levels varying
% with time

% *Incidentally, this is a great system for which to use simbiology

% Time here is in minutes 

plot(prottimes,protexps);
title('Protein Expression');
xlabel('Time (min)');
ylabel('Expression level');
shg;


%%
% What appears to be happening here becomes even clearer when we identify
% which protein is which
legend(prot_labels{:}, 'location', 'NorthEastOutside');
shg;


%% Signal 3: EEG Data
% The EEG data we will cover in more detail, but as general info:

% What is an EEG? Electroencephalography - as Todd and Justin covered - 
% records voltage changes at the scalp surface so we can infer about underlying
% activity in the brain.

% High temporal but low spatial resolution (good for epilepsy, sleep)
% We typically focus on neural oscillations for EEG.

% One helpful analogy is think of the scalp electrodes as buoys floating on
% water in a lake. 

% There are some people throwing stones on the other side of the lake.

% The rise and fall of the buoys reflects the waves caused by the stones.
% Let's infer where and when the stones dropped from the buoy movement.

% Analogously let us infer about oscillatory activity in the neural
% circuits that cause changes observed on EEG.

% EEG01 and EEG02 are each one minute EEG recordings from an epileptic pediatric patient
% The kid had a seizure during one and not the other minute

% Data courtesy of PhysioNet

% Plot 23 channel EEG01
figure;
subplot(2,1,1);
plot(EEGtimes,EEG01);
title('EEG01');
xlabel('Time (sec)');
ylabel('Voltage (uV)');
legend(ch_labels{:}, 'location', 'NorthEastOutside');
shg;

% Plot 4 channel EEG01 
subplot(2,1,2);
plot(EEGtimes,EEG01(1:4,:));
title('EEG01');
xlabel('Time (sec)');
ylabel('Voltage (uV)');
legend(ch_labels{1:4}, 'location', 'SouthEastOutside');
shg;


%% % Signal Processing

% What is a signal? 
% Any function that "conveys information about the behavior or attributes of some phenomenon"
% Particularly interested in time-variant signals like the ones we
% discussed, but an image is also a signal that varies over space rather
% than time.

%% Power of a signal
% The energy of a signal is equal to sum of squares of samples
% The power of a signal is equal to sum of squares of samples divided by signal length (means squared)
energ = sum(temps.^2);
pEnerg = energ/length(temps)
pRMS = rms(temps)^2

% pwelch returns power spectral density estimate, divides up by frequency
% and gives psd by using bandpower()
figure;
pwelch(temps,[],[],[],365*2);
xlabel('Cycle/Year');

%% Signal cross-correlation

% Our goal is often to identify significant relationships between signals. A useful way to do that
% can be to compute and display cross-correlations between channels.


% Cross-correlation measures the relatedness of two signals by trying to
% shift one so that it perfectly aligns with the other.
% (Some of you may know about convolution, a related method for 2 continuous
% functions, which takes the integral of the functions multiplied with one reversed
% and shifted).

% In the protein example, let us find possible relationships between these signals.

% Creating a correlation matrix
[cc,p] = corrcoef(protexps');

% Let us list which signal pairs are significantly correlated (at alpha =
% 0.05), removing redundant pairs
[ind1,ind2] = find(0 < p - tril(p) & p - tril(p) < 0.05);
sigpairs = [prot_labels(ind1); prot_labels(ind2)]

% Visualizing the correlation matrix
imagesc(cc);
colorbar;
title('Protein Expression Cross-Correlations');
set(gca,'XLim',[0.5 4.5],'XTick',1:4,'YLim',[0.5 4.5],'YTick',1:4,...
    'XTickLabel',prot_labels,'YTickLabel',prot_labels);


% Note that although cross-correlation is conveniently used to detect
% underlying network dynamics in neurophysiology, it tells us nothing about
% causality!!

%% Spectral Analysis
% Fourier said give me any function and I can represent it as a sum of
% composite frequencies.

% FFT is a fast algorithm to compute the FT, which translates time domain
% to frequency domain.

% When you actually run the transform, the absolute value is how much is
% present in the original function. The complex argument is phase offset.

% Let us FFT the temperature data since that seemed clearly oscillatory
tempfft = fftshift(fft(temps));

% In the complex components of tempfft, there is an embedded magnitude and phase
tempfftmag = abs(tempfft);
tempfftphase = angle(tempfft);

% This is our plot of signal energy (2-sided)
figure;
plot(tempfftmag.^2);
title('Signal Energy')
xlabel('Cycles/Year')
shg;

% This is our plot of signal power (1-sided)
figure;
pwelch(temps,[],[],[],365*2);
xlabel('Cycles/Year')
shg;

% A spectrogram of the signal
figure;
spectrogram(temps,[],[],[],365*2);
xlabel('Cycles/Year')

%% EEG Frequency Analysis

% If we want to identify interesting events on EEG, we can look at signal and 
% spectrogram at the same time.
simuleegspec(EEGtimes,EEG01(1,:),Fs);



%%  In-Class Exercise 
% Visualize the two EEGs. We will try to figure out which one contains the
% seizure.
% Cross correlation during seizure versus during normal activity:
% Find the  channel cross-correlations in each set that are significant at a =
% 0.05.
% What does a significant value mean?
% Visualize the correlation matrix between channels
% What's the difference between the two?

% FFT the two datasets and compare the spectrograms.
% What does this signify?

% Which dataset do you think is the one with the seizure?

% Can you visualize the cross-correlogram for key parts of the seizure event?

figure;
subplot(2,1,1);
plot(EEGtimes,EEG01);
title('EEG01');
xlabel('Time (sec)');
ylabel('Voltage (uV)');
legend(ch_labels, 'location', 'NorthEastOutside');
shg; 
subplot(2,1,2);
plot(EEGtimes,EEG02);
title('EEG02');
xlabel('Time (sec)');
ylabel('Voltage (uV)');
legend(ch_labels, 'location', 'SouthEastOutside');
shg;

simuleegspec(EEGtimes,EEG01(1,:),Fs);
simuleegspec(EEGtimes,EEG02(1,:),Fs);

figure;
subplot(2,1,1)
imagesc(corrcoef(EEG01'));
colorbar;
title('EEG01 Cross-Correlations');
set(gca,'XLim',[0.5 23.5],'XTick',1:23,'YLim',[0.5 23.5],'YTick',1:23,...
    'YTickLabel',ch_labels);


subplot(2,1,2)
imagesc(corrcoef(EEG02'));
colorbar;
title('EEG02 Cross-Correlations');
set(gca,'XLim',[0.5 23.5],'XTick',1:23,'YLim',[0.5 23.5],'YTick',1:23,...
    'YTickLabel',ch_labels);

seizure = EEG02(:,Fs*35:Fs*50);
figure;
imagesc(corrcoef(seizure'));
colorbar;
title('Seizure Cross-Correlations');
set(gca,'XLim',[0.5 23.5],'XTick',1:23,'YLim',[0.5 23.5],'YTick',1:23,...
    'YTickLabel',ch_labels);

% EEG01 is normal EEG.
% EEG02 contains the seizure!

%% %% Extra Topics

%% % Aliasing in temporal and image signals 

% Aliasing refers to distortion in the signal that results from sampling
% limitations. This can be in signals that are either temporal or spatial sampled.

%% Aliasing Temperature Data

% What if we downsampled our temperature data by a factor of 2 (taking
% every other data point)?
tempsdown = temps(2:2:length(temps));
temptimesdown = temptimes(2:2:length(temps));

figure;
subplot(2,1,1);
plot (temptimesdown, tempsdown);
title ('Temperature in Boston');
xlabel('Time (days)');
ylabel('Temp (F)');
shg;

subplot(2,1,2);
plot (temptimes, temps);
title ('Temperature in Boston');
xlabel('Time (days)');
ylabel('Temp (F)');
shg;

figure;
subplot(2,1,1);
spectrogram(tempsdown,[],[],[],365);
xlabel('Cycles/Year')
shg;

subplot(2,1,2);
spectrogram(temps,[],[],[],365*2);
xlabel('Cycles/Year')
shg;

% Note how we lose information on the daily oscillation if we only sample every
% other timepoint and could not differentiate.

% The fact that the highest frequency we can test for is 365 cycles/year is
% based on our sampling frequency of 365*2 samples/year.

% The frequency of interest (365 cycles/year) is the Nyquist frequency.
% The sampling rate for that frequency (2x) is the Nyquist rate.


%% Aliasing An Image

% Image example
% Creating pic 1 = original
phelps1 = imread('Phelps1.png');
phelps1 = rgb2gray(phelps1);

lochte1 = imread('Lochte1.png');
lochte1 = rgb2gray(lochte1);
phelpsclean = imresize(phelps1, size(lochte1));
figure;
imshow([lochte1 phelpsclean])

% Adding Lochte (noise) to winning signal Phelps
indices = 4:4:size(phelpsclean,2);
noise = lochte1;
phelpsnoisy = phelpsclean;
phelpsnoisy(:,indices) = noise(:,indices);

% Here we create downsampled versions of both pics (taking every other
% datapoint, you can change 2 to modify downsampling factor) 
phelpscleandown = phelpsclean(4:4:end,4:4:end);
phelpsnoisydown = phelpsnoisy(4:4:end,4:4:end);
lochte1down = lochte1(4:4:end,4:4:end);

% Here are the original pics, note that they appear different 
% due to the bars  
figure;
imshow([phelpsclean,phelpsnoisy,lochte1]);
title('original clean Phelps      |    original noisy Phelps |  original Lochte');

% Here are the downsampled pics, now they appear the same
figure;
imshow([phelpscleandown,phelpsnoisydown,lochte1down])
title('downsampled clean Phelps   |   downsampled noisy Phelps |  downsampled Lochte');

% The loss of information such that we cannot distinguish the two pictures
% is aliasing. We generally do not want to have this happen (assuming the
% difference between the pictures is important to us).

% Fortunately for us in physiology, we tend to consider the higher
% frequencies noise so aliasing is less of a problem at <kHz sampling (compared to
% say audio or image processors).

