%Cem Rýfký Aydýn 2013800054
%CmpE58P - MACHINE LISTENING: Homework II
%Question 1

%All past data get cleared
clear all;
close all;
clc;

%All the vowels are stored in the following directory called 'vowels'
dirName = 'vowels\';
%If one wants to see the subplots for the data-set "vowels_Web", s/he has
%to uncomment the following
%dirName = 'vowels_Web\';
files = dir(dirName);

fileN = [];

%We iterate over the files in the directory through the below loop
for file = files';
    
    if strcmp(file.name, '.') == 0 && strcmp(file.name, '..') == 0
        fileN = [fileN; char(strcat(strcat(dirName, char(file.name))))];
    end
    
end

allData = cellstr(fileN);

for u = 1:length(allData)
    
    fileMv = char(allData(u));
    
    
    %%
    
    %__PLOTS_______________________________________________________________
    
    [data,freks_,bits] = wavread(fileMv);
    time_ = 0:1/freks_:(length(data)-1)/freks_;
    figure(u);
    subplot(2, 3, 1)
    plot(time_, data);
    title(strcat('Waveform for "', fileMv(8:end), '"')); 
    xlabel('Time (secs)');
    ylabel('Amplitude');
    
    %Spectrogram being plotted
    N_ = 512;
    subplot(2,3,3);
    rounded =round(N_ * .25 * 0.91);
    spectrogram(data, hamming(N_/4), rounded, N_, freks_);
    title('Spectrogram');
    xlabel('Frequency (in Hz)');
    ylabel('Time (secs)');
    
    subplot(2,3,5)


    %Fourier transformation
    Four_ = fft(hamming(length(data)) .* data);

    highFreq = 2500 * 2 * length(Four_)/freks_;
    f = (0:highFreq) * freks_ / length(Four_);
    subplot(2, 3, 5);
    
    % The built-in function "cceps" could also have been used
    % Cepstrum is the Fourier transform of the log spectral density
    C=fft(log(abs(Four_) + 0.00000001));
    
    
    xRange1 = 1;
    xRange2 = length(data);
    quef = (xRange1:xRange2)/freks_;
    plot(quef, abs(C(xRange1:xRange2)));
    title('Cepstrum values');
    xlabel('Quefrency (secs)');
    ylabel('Amplitude');
    
    %A short interval as an example
    x1_ = 0.08;
    x2_ = 0.20;
    subplot(2,3,2)
    plot(time_,data)
    xlim([x1_ x2_])
    title('Waveform of time')
    xlabel('Time (secs)')
    ylabel('Amplitude')
    
               
    xRange2 = freks_ / 50; %50Hz
    xRange1 = freks_ / 1000;  %1000Hz
    
    Four_ = fft(hamming(length(data)) .* data);

    freqCount = 5000;
    highFreq = length(Four_) / freks_ * freqCount;
    frequ_ = (0:highFreq) * freks_/length(Four_);
    
    subplot(2, 3, 2);
    plot(frequ_, log10(abs(Four_(1:length(frequ_)))+ 0.0000001) * 20);
    title('Spectrum values');
    xlabel('Frequency (in Hz)');
    ylabel('Magnitude (in dB)');
    
    %Extracting 30 ms of the waveform

    subplot(2, 3, 4)
    xstart = 0.050;
    plot(time_,data)
    xlim([xstart (0.0300 + xstart)])

    samples=data(xstart*freks_:xstart * freks_ + .03 * freks_);
    
    
    
    % cepstrum being calculated
    Ceps_=fft(log(abs(Four_)+eps));
    % cepstrum values captured in a short interval are being plotted again,
    % so as to compare it with the overall cepstral values in the whole
    % time domain
    quef = (xRange1:xRange2) ./ freks_;
    subplot(2, 3, 4);
    
    plot(quef, abs(Ceps_(xRange1:xRange2)));
    title('Cepstrum values');
    xlabel('Quefrency (secs)');
    ylabel('Amplitude');

    
     %Spectral envelope and magnitude
    frLen = 512;
    frCnt = 8000;
    freq = 0:1/N_ * frCnt:frCnt - 1;
    magnit_ = abs(fft(samples, frLen * 2));
    magnit_ = magnit_(1:frLen);

    subplot(2, 3, 6)
    plot(freq,db(magnit_/(frLen * 2)));      % divide by 1024 to normalize after fourier transform
    grid on; axis tight; title('');
 
    hold on;
    
    
    %Smoothing envelope being drawn below
    p = 0.001 * freks_;
    p = p + 4;
    [lpcA, lpcG] = lpc(samples, p);
    lspec = freqz(lpcG, lpcA, freq, freks_);
    subplot(2,3,6)
    plot(freq, log10(abs(lspec)) * 20, 'r');
    xlabel('Frequency (in Hz)'); ylabel('||X|| (in dB)');
    title('Spectral envelope: red, magnitude: blue');     axis tight;

    
    
    %%
    
    %FORMANT FREQUENCY SECTION
    clear y Fs;
    
    
    [data, Fs] = wavread(fileMv);
    
    sound(data, Fs);
    
    %data = data(round(length(data) / 12):round(length(data) - length(data) / 12));
    NFFT = 128;
    
    len_ = 105;
    
    overlappNo = 95;
    
    
    
    
    
    xMod = data .* hamming(length(data));
    
    %Filtering is implemented below
    xMod = filter(1, [1 0.63], xMod);
    
    %How many formant frequencies are to extract is stated through the
    %below statement
    formantNumb = 2;
    %LPC method being implemented below
    outp_ = lpc(xMod, formantNumb * 2 + 2);
    roots_ = roots(outp_);
    
    %Mathematical equations being performed for calculating the formant
    %frequencies
    roots_ = roots_(imag(roots_)>=0);
    angulars_ = atan2(imag(roots_),real(roots_));
    
    [freks,ind_] = sort((Fs/2/pi) .* angulars_);
    bandwidth = -.5 * log(abs(roots_(ind_))) *(Fs/2/pi);
    
    formantVals = [];
    formNo = 1;

    it_ = 1;
    if strcmp( dirName, 'vowels_Web\')
        cmpBand_ = 10000;
    else
        cmpBand_ = 400;
    end
    while it_ <= length(freks)
        if bandwidth(it_) < cmpBand_
            if (freks(it_) > 90)
                formantVals(formNo) = freks(it_);
                formNo = formNo+1;
            end
        end
        it_ = it_ + 1;
    end
    fileMv(8:end)
    formantVals(end:-1:end-1)
end




