% get a section of vowel
[x,fs]=wavread('i.wav');

clc;
close all;

signal = x;

Fs = fs;

clc;
close all;

%Normal FFT
y=signal;
figure();
N=2*2048;T=N/Fs;
sig_f=abs(fft(y(1:N)',N));
sig_n=sig_f/(norm(sig_f));
freq_s=(0:N-1)/T;
plot(freq_s(2:250),sig_n(2:250));title('FFT of Original Signal');


%Envelope Detection based on Low pass filter and then FFT
[a,b]=butter(2,0.1);%butterworth Filter of 2 poles and Wn=0.1
%sig_abs=abs(signal); % Can be used instead of squaring, then filtering and
%then taking square root
sig_sq=2*signal.*signal;% squaring for rectifing 
%gain of 2 for maintianing the same energy in the output
y_sq = filter(a,b,sig_sq); %applying LPF
y=sqrt(y_sq);%taking Square root
%advantages of taking square and then Square root rather than abs, brings
%out some hidden information more efficiently
figure();
N=2*2048;T=N/Fs;
sig_f=abs(fft(y(1:N)',N));
sig_n=sig_f/(norm(sig_f));
freq_s=(0:N-1)/T;
plot(freq_s(2:250),sig_n(2:250));title('Envelope Detection: LPF Method');



%Envelope Detection based on Hilbert Transform and then FFT
analy=hilbert(signal);
y=abs(analy);
figure();
N=2*2048;T=N/Fs;
sig_f=abs(fft(y(1:N)',N));
sig_n=sig_f/(norm(sig_f));
freq_s=(0:N-1)/T;
plot(freq_s(2:250),sig_n(2:250));title('Envelope Detection : Hilbert Transform')

