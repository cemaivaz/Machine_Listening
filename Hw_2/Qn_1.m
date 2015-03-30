close all;
clear all;
clc;
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


figure(2);
% Read the data back into MATLAB, and listen to audio.

[y, Fs, nbits, readinfo] = wavread('i.wav');

sound(y, Fs);

%Below, spectral density is implemented

len = 100;

noverlap = 90;

NFFT = 128;

[y_,freq,time,p] = spectrogram(y,len,noverlap,NFFT,Fs);

surf(time,freq,10*log10(abs(p)),'EdgeColor','none');

axis xy; axis tight; colormap(jet); view(0,90);

xlabel('Time');

ylabel('Frequency (in terms of Hz)');

tSpl = 1/Fs;

stInd = round(0.2/tSpl)

endInd = round(0.8/tSpl)

%vals = y(stInd:endInd);

vals = y(1:size(y, 1));

x1 = vals.*hamming(length(vals));

preemph = [1/4 1/4 1/4 1/4];

x1 = filter(1,preemph,x1);


formNo = 3;
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

formantVals

%%

% Create WAV file in the current folder.

load handel.mat

figure(3);

fileSnd = 'handel.wma';

wavwrite(y, Fs, fileSnd)

%clear y Fs

% Read the data back into MATLAB, and listen to audio.

[x,rate]=wavread('i.wav');

x=x(1:500);

x=x.*hamming(500); %create hamming windows

p=12;

a=lpc(x,p); %implement lpc method

nfft=512;

freq=linspace(0,rate/2,nfft/2+1);

FFT_log = 20 * log10 ( abs( fft( x, nfft ) ) ); %decibel is measured 

%terms of logarithm

FFT_log = FFT_log(1:nfft/2+1); %FFT implemented for the above audio 

plot(freq,FFT_log,'b');

LPC_log = -20 * log10 ( abs( fft( a, nfft ) ) );

LPC_log = LPC_log(1:nfft/2+1); %superimposed LPC spectrum

hold on;

plot(freq,LPC_log-3,'r');

xlabel('Frequency (in terms of Hertz)')

ylabel('Magnitude (in terms of Decibel)')

hold off