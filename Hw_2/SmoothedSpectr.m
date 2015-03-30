%Time Waveform
clear all;
close all;
[x,fs,bits] = wavread('i.wav');
t = 0:1/fs:(length(x)-1)/fs;
figure(1);
subplot(3,2,1)
plot(t,x);
title('Time Waveform of the Word example'); xlabel('Time(s)');ylabel('Amplitude')  


%Narrowband Spectrogram and WideBand Spectrogram
N = 512;
subplot(3,2,3)
spectrogram(x,hamming(N/4),round(0.9*N/4),N,fs);
title('NarrowBand Spectrogram');

N = 512;
subplot(3,2,5)
win = hamming(N);
spectrogram(x,win,round(0.97*N),N,fs);
title('Wideband spectrogram');

%--------------------------------------------------------------------------

%Time Waveform of Phoneme 
xa = 0.06;                 % These values where varied according to the respective waveform of the word being plotted
xb = 0.17;
subplot(3,2,2)
plot(t,x)
xlim([xa xb])
title('Time Waveform of Phoneme ')
xlabel('Time (s)')
ylabel('Amplitude')

%Extracting 30 ms of the waveform 
xstart = 0.05;
subplot(3,2,4)
plot(t,x)
xlim([(xstart) (xstart+.030)])
title('Time Waveform of 30 ms of the Phoneme')
xlabel('Time (s)')
ylabel('Amplitude')
sample=x(xstart*fs:(xstart+.030)*fs);

%Magnitude Spectrum and Linear Prediction Spectral Envelope
mag = abs(fft(sample,1024));
mag = mag(1:512);     
freq = 0:8000/N:8000-1;
subplot(3,2,6)
plot(freq,db(mag/1024));      % divide by 1024 to normalize after fourier transform
axis tight; grid on; title('Magnitude Spectrum and Smoothed Spectral Envelope');
xlabel('Frequency (Hz)'); ylabel('|X(f)| (dB)');
hold on


%Smoothing - Linear Prediction
p=fs/1000 + 4;
[a,g]=lpc(sample,p);
lspec = freqz(g,a,freq,fs);
subplot(3,2,6)
plot(freq, 20*log10(abs(lspec)),'k');
axis tight;
title('Magnitude Spectrum in Blue & Smoothed Spectral Envelope in Black');

%3D PLOT as a function of Time, Frequency and Power Density
figure(2)
fD = 0:8000/N:8000-1;
fD = fD';
[S,F,T,W] = spectrogram(x,hamming(N),floor(N/9),fD,fs);
[a,b]=meshgrid(-2:1:2,-2:1:2);
surf(freq,T,db(abs(W'))); hold on;
title('3-D plot of recorded word')
xlabel('Frequency,  Hz');
ylabel('Time, s');
zlabel('Power Spectral Density, dB');
h = surf(freq,T,db(abs(W')));
set(h,'edgecolor','none')