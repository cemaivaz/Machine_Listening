clear
clc

% Read the data back into MATLAB, and listen to audio.

fileSnd = '03_CelloTaksimi_pt1.wav';

[y, Fs, nbits, readinfo] = wavread(fileSnd);
%If it had been native, unnormalized data would have been taken into
%account, with the values thereof ranging from -32,000 to +32,000 for the
%bit count of 16

%sound(y, Fs);

%Below, spectral density is implemented

%How many windows?
len = 100;

noverlap = 50;

%The below constant defines the length (# of samples) of each window
NFFT = 128;

[y_,freq,time,p] = spectrogram(y,len,noverlap,NFFT,Fs);


surf(time,freq,10*log10(abs(p)),'EdgeColor','none');

axis xy; axis tight; colormap(jet); view(0,90);

xlabel('Time');

ylabel('Frequency (in terms of Hz)');

windLength = Fs * 0.04;
windLength = 64;

% onsetdet(f,win_length,thre,range,multi,shift,showplot,realtime)
[pos, v] = onsetdet(y,windLength,3,10,2,12,1,length(y) / Fs);
%[pos, v] = onsetdet(y,windLength,3,10,2,12,1,length(y) / Fs);
%Sonuncu parametre ile onset'lerin kaçýncý saniyede oluþtuðunu görüyorsun