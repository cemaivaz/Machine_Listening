clear all;
close all;
clc;
files = dir('vowels');

fileN = [];

for file = files';
    
    char(file.name)
    if strcmp(file.name, '.') == 0 && strcmp(file.name, '..') == 0
        fileN = [fileN; char(strcat(strcat('vowels\', char(file.name))))];
    end
    
end

allData = cellstr(fileN);

for u = 1:length(allData)
    
    
    fileMv = allData(u);
    fileMv = char(fileMv);

end



figure(8);
% get a section of vowel
[x,fs]=wavread('i.wav');
ms1=fs/1000;                 % maximum speech Fx at 1000Hz
ms20=fs/50;                  % minimum speech Fx at 50Hz

% ms1=1;
% ms20=length(x);
%
% plot waveform
t=(0:length(x)-1)/fs;        % times of sampling instants
subplot(3,1,1);
plot(t,x);
legend('Waveform');
xlabel('Time (s)');
ylabel('Amplitude');
%
% do fourier transform of windowed signal
Y=fft(x.*hamming(length(x)));
%
% plot spectrum of bottom 5000Hz
hz5000=5000*length(Y)/fs;
f=(0:hz5000)*fs/length(Y);
subplot(3,1,2);
plot(f,20*log10(abs(Y(1:length(f)))+eps));
legend('Spectrum');
xlabel('Frequency (Hz)');
ylabel('Magnitude (dB)');
%
% cepstrum is DFT of log spectrum
C=fft(log(abs(Y)+eps));
%
% plot between 1ms (=1000Hz) and 20ms (=50Hz)
q=(ms1:ms20)/fs;
subplot(3,1,3);
plot(q,abs(C(ms1:ms20)));
legend('Cepstrum');
xlabel('Quefrency (s)');
ylabel('Amplitude');

% figure(8);
% c = cceps(s2);
%
% plot(t,c)
% xlabel('Time (s)')
% title('Complex cepstrum')

%%

%SMOOTHED PLOT

%Time Waveform

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

% do fourier transform of windowed signal
Y=fft(x.*hamming(length(x)));
%
% plot spectrum of bottom 5000Hz
hz5000=5000*length(Y)/fs;
f=(0:hz5000)*fs/length(Y);
subplot(3,2,3);

%
% cepstrum is DFT of log spectrum
C=fft(log(abs(Y)+eps));
%
% plot between 1ms (=1000Hz) and 20ms (=50Hz)
ms1 = 1;
ms20 = length(x);
q=(ms1:ms20)/fs;
plot(q,abs(C(ms1:ms20)));
title('Cepstrum');
xlabel('Quefrency (s)');
ylabel('Amplitude');
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
plot(freq, 20*log10(abs(lspec)),'r');
axis tight;
title('Magnitude Spectrum in Blue & Smoothed Spectral Envelope in Black');




%%

%FORMANT FREQUENCIES


% Read the data back into MATLAB, and listen to audio.

[y, Fs, nbits, readinfo] = wavread('i.wav');

sound(y, Fs);

%Below, spectral density is implemented

len = 100;

noverlap = 90;

NFFT = 128;

[y_,freq,time,p] = spectrogram(y,len,noverlap,NFFT,Fs);



tSpl = 1/Fs;

stInd = round(0.2/tSpl)

endInd = round(0.8/tSpl)

%vals = y(stInd:endInd);

vals = y;

x1 = vals.*hamming(length(vals));

preemph = [1/4 1/4 1/4 1/4];

x1 = filter(1,preemph,x1);


formNo = 2;
A = lpc(x1, formNo * 2 + 2);


roots = roots(A);

roots = roots(imag(roots)>=0);

angz = atan2(imag(roots),real(roots));

[freqs,index] = sort(angz.*(Fs/(2*pi)));

bw = -1/2*(Fs/(2*pi))*log(abs(roots(index)));

i = 1;

[c3, ind] = max(size(freqs));

for j = 1:c3
    
    if (freqs(j) > 90 && bw(j) <400)
        
        formantVals(i) = freqs(j);
        
        i = i+1;
        
    end
    
end

formantVals(end:-1:end - 1)