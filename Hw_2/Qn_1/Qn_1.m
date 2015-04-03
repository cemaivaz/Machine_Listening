%Cem Rýfký Aydýn 2013800054
%CmpE58P - MACHINE LISTENING: Homework II
%Question 3

%All past data get cleared
clear all;
close all;
clc;

%All the vowels are stored in the following directory called 'vowels'
files = dir('vowels');

fileN = [];

%We iterate over the files in the directory through the below loop
for file = files';
    
    if strcmp(file.name, '.') == 0 && strcmp(file.name, '..') == 0
        fileN = [fileN; char(strcat(strcat('vowels\', char(file.name))))];
    end
    
end

allData = cellstr(fileN);

for u = 1:length(allData)
    
    
    fileMv = char(allData(u));
    
    
    %%
    
    %Graphics section
    
    %Time Waveform
    
    [data,fs] = wavread(fileMv);
    t = 0:1/fs:(length(data)-1)/fs;
    figure(u);
    subplot(2,3,1)
    plot(t,data);
    title(strcat('The waveform - "', fileMv(8:end), '"')); 
    xlabel('Time (secs)');
    ylabel('Amplitude')
    
    
    %Narrowband Spectrogram and WideBand Spectrogram
    windLen = 512;
    subplot(2,3,3)
    rounded = round(1 / 4 * 0.89 * windLen)
    spectrogram(data,hamming(windLen * 1 / 4),rounded,windLen,fs);
    title('Spectrogram');
    
    % Fourier transform is implemented below
    Y=fft(hamming(length(data)) .* data);
   
    %Plotting the cepstrum values
    subplot(2,3,5);
    
    % The built-in function "cceps" could have also been used here
    % Cepstrum can be defined as the fft taking the log-spectrum as input
    C=fft(log(eps + abs(Y)));
    
    
    firstVal = 1;
    secVal = length(data);
    q=(firstVal:secVal)/fs;
    plot(q,abs(C));
    title('Cepstrum values');
    xlabel('Quefrency (sec)');
    ylabel('Amplitude');

    
    %Time Waveform of Phoneme
    xa = 0.06;                 % These values where varied according to the respective waveform of the word being plotted
    xb = 0.17;
    subplot(2,3,2)
    plot(t,data)
    xlim([xa xb])
    title('Time Waveform of Phoneme ')
    xlabel('Time (s)')
    ylabel('Amplitude')
    
    firstVal=fs/1000;                 % maximum speech Fx at 1000Hz
    secVal=fs/50;                  % minimum speech Fx at 50Hz
    
    Y=fft(data.*hamming(length(data)));
    %
    % plot spectrum of bottom 5000Hz
    hz5000=5000*length(Y)/fs;
    f=(0:hz5000)*fs/length(Y);
    subplot(2,3,2);
    plot(f,20*log10(abs(Y(1:length(f)))+eps));
    title('Spectrum');
    xlabel('Frequency (Hz)');
    ylabel('Magnitude (dB)');
    
    %Extracting 30 ms of the waveform
    xstart = 0.05;
    subplot(2,3,4)
    plot(t,data)
    xlim([(xstart) (xstart+.030)])
    title('Time Waveform of 30 ms of the Phoneme')
    xlabel('Time (s)')
    ylabel('Amplitude')
    sample=data(xstart*fs:(xstart+.030)*fs);
    
    
    
    % cepstrum is DFT of log spectrum
    C=fft(log(abs(Y)+eps));
    %
    % plot between 1ms (=1000Hz) and 20ms (=50Hz)
    q=(firstVal:secVal)/fs;
    subplot(2,3,4);
    plot(q,abs(C(firstVal:secVal)));
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
    freq = 0:8000/windLen:8000-1;
    subplot(2,3,6)
    plot(freq,db(mag/1024));      % divide by 1024 to normalize after fourier transform
    axis tight; grid on; title('Magnitude Spectrum and Smoothed Spectral Envelope');
    xlabel('Frequency (Hz)'); ylabel('|X(f)| (dB)');
    hold on
    
    
    %Smoothing - Linear Prediction
    p=fs/1000 + 4;
    [a,g]=lpc(sample,p);
    lspec = freqz(g,a,freq,fs);
    subplot(2,3,6)
    plot(freq, 20*log10(abs(lspec)),'r');
    axis tight;
    title('Magnitude Spectrum in Blue & Smoothed Spectral Envelope in Black');
    
    
    
    
    %%
    
    %FORMANT FREQUENCY SECTION
    clear y Fs;
    
    
    [data, Fs] = wavread(fileMv);
    
    %The file is listened through the below code
    sound(data, Fs);
    
    %The features of frames are given below
    overlappNo = 95;
    
    len_ = 110;
    
    NFFT = 128;
    
    %Hamming window is utilized
    xMod = data.*hamming(size(data, 1));
    
    %Filtering is performed below
    xMod = filter(1,[1 0.7],xMod);
    
    %In this homework, we were supposed to take the number of formants as 2
    noFormants = 2;
    
    %LPC method
    output = lpc(xMod,noFormants *3 + 2);
    
    %In order to extract the formant frequencies, the mathematical
    %equations below are performed
    roots_ = roots(output);
    
    
    roots_ = roots_(imag(roots_)>=0);
    angulars_ = atan2(imag(roots_),real(roots_));
    
    [freks,ind_] = sort(angulars_.*(Fs/2/pi));
    bandwidth_ = -1/2*log(abs(roots_(ind_)))*(Fs/(2*pi));
    
    formantVals = [];
    formNo = 1;

    it_ = 1;
    while it_ <= length(freks)
        if bandwidth_(it_) < 400
            if (freks(it_) > 90)
                formantVals(formNo) = freks(it_);
                formNo = formNo+1;
            end
        end
        it_ = it_ + 1;
    end
    fileMv(8:end)
    %The top 2 formant frequencies
    formantVals(end:-1:end-1)
end




