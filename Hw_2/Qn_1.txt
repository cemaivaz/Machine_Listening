clear all;
close all;
clc;
files = dir('vowels');

fileN = [];

for file = files';
    
    if strcmp(file.name, '.') == 0 && strcmp(file.name, '..') == 0
        fileN = [fileN; char(strcat(strcat('vowels\', char(file.name))))];
    end
    
end

allData = cellstr(fileN);

for u = 1:length(allData)
    
    
    fileMv = allData(u);
    fileMv = char(fileMv);
    
    %%
    
    %SMOOTHED PLOT
    
    %Time Waveform
    
    [x,fs,bits] = wavread(fileMv);
    t = 0:1/fs:(length(x)-1)/fs;
    figure(u);
    subplot(3,2,1)
    plot(t,x);
    title(strcat('Time Waveform of "', fileMv(8:end), '"')); xlabel('Time(s)');ylabel('Amplitude')
    
    
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
    subplot(3,2,5);
    
    % cceps could have also been used
    % cepstrum is DFT of log spectrum
    C=fft(log(abs(Y)+eps));
    %5
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
    
    ms1=fs/1000;                 % maximum speech Fx at 1000Hz
    ms20=fs/50;                  % minimum speech Fx at 50Hz
    
    Y=fft(x.*hamming(length(x)));
    %
    % plot spectrum of bottom 5000Hz
    hz5000=5000*length(Y)/fs;
    f=(0:hz5000)*fs/length(Y);
    subplot(3,2,2);
    plot(f,20*log10(abs(Y(1:length(f)))+eps));
    title('Spectrum');
    xlabel('Frequency (Hz)');
    ylabel('Magnitude (dB)');
    
    %Extracting 30 ms of the waveform
    xstart = 0.05;
    subplot(3,2,4)
    plot(t,x)
    xlim([(xstart) (xstart+.030)])
    title('Time Waveform of 30 ms of the Phoneme')
    xlabel('Time (s)')
    ylabel('Amplitude')
    sample=x(xstart*fs:(xstart+.030)*fs);
    
    
    
    % cepstrum is DFT of log spectrum
    C=fft(log(abs(Y)+eps));
    %
    % plot between 1ms (=1000Hz) and 20ms (=50Hz)
    q=(ms1:ms20)/fs;
    subplot(3,2,4);
    plot(q,abs(C(ms1:ms20)));
    title('Cepstrum');
    xlabel('Quefrency (s)');
    ylabel('Amplitude');
    
    
    %
    % do fourier transform of windowed signal
    
    %
    
    %
    
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
    clear y Fs;
    
    % Read the data back into MATLAB, and listen to audio.
    
    [y, Fs, nbits, readinfo] = wavread(fileMv);
    
    sound(y, Fs);
    
    %Below, spectral density is implemented
    
    len = 100;
    
    noverlap = 90;
    
    NFFT = 128;
    
    [y_,freq,time,p] = spectrogram(y,len,noverlap,NFFT,Fs);
    
    x = y;
    
    x1 = x.*hamming(length(x));
    
    preemph = [1 0.63];
    x1 = filter(1,preemph,x1);
    
    
    A = lpc(x1,8);
    rts = roots(A);
    
    
    rts = rts(imag(rts)>=0);
    angz = atan2(imag(rts),real(rts));
    
    [frqs,indices] = sort(angz.*(Fs/(2*pi)));
    bw = -1/2*(Fs/(2*pi))*log(abs(rts(indices)));
    
    formants = [];
    nn = 1;
    for kk = 1:length(frqs)
        if (frqs(kk) > 90 && bw(kk) <400)
            formants(nn) = frqs(kk);
            nn = nn+1;
        end
    end
    fileMv(8:end)
    formants
end




