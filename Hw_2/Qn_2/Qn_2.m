%Cem R�fk� Ayd�n 2013800054
%CmpE58P - MACHINE LISTENING: Homework II
%Question 2

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
    
    
    p = [];
    [data fs] = wavread(fileMv);
    
    timeData = length(data) / fs;
    
    
    %My implementation for YIN algorithm is given below. I utilized
    %interpolating parabolic equation for it.
    
    dataSize = length(data);
    
    thre_ = 0.11; %yin threshold
    range_ = 256;
    
    windLen = range_;
    skip = windLen; %It denotes the frame period
    
    
    frameNo = floor((dataSize - range_ - windLen) / skip); %The number of frames
    
    
    period_ = ones(1, frameNo);
    p_ = ones(1, frameNo); %pwr
    
    
    data = convmtx(data, range_ + 1.); % Data gets shifted
    
    k = 1;
    data = data(range_:end-range_,:);

    while k <= frameNo
        
        %The frame offset is defined below
        st_ = k * skip - skip; 
        xSub = data(st_ + 1:st_ + windLen,:);
        repAmount = range_ + 1;
        %The difference is stated below (mean square)
        dist_ = .5 * mean(power(xSub - repmat(xSub(:,1),1,repAmount), 2)); 
        %cumM_ = cumulative mean
        cumM_ = dist_(2:length(dist_)) ./ (cumsum(dist_(2:length(dist_))) ./ (1:(repAmount - 1)));
        
        % Interpolation getting implemented below
        positionsMin=1:length(cumM_);
        
        firstX = cumM_(1:length(cumM_)-2);
        secX = cumM_(2:length(cumM_)-1);
        thirdX = cumM_(3:length(cumM_));
        %Parabolic coefficients in accordance with the algorithm I developed
        %are given below
        a_ = .5 * (firstX + thirdX - 2. * secX);
        b_ = .5 * (thirdX - firstX);
        sh_ = -b_ ./ a_ * .5;
        output = secX - b_ .^2 ./ a_ .* 0.25;
        
        
        ind_ = find(secX<firstX & secX<thirdX) + 1.;
        cumM_(ind_) = output(-1. + ind_);
        positionsMin(ind_) = positionsMin(-1. + ind_) + sh_(-1. + ind_);
        
        %Minimum values are getting extracted below
        
        a_ = zeros(length(cumM_));
        for o_ = 1:length(cumM_)
            a_(o_) = (cumM_(o_) < thre_);
        end
        
        %Global minimum or local minimum are getting taken into account
        %below
        if length(find(a_)) == 0
            
            [minVal, indMin] = min(cumM_);
        else
            bTmp = find(a_);
            bTmp = bTmp(1);
            b_ = min(bTmp);
            c_ = min(b_*2., length(a_));
            [minVal, indMin] = min(cumM_(b_:(c_-1)));
            indMin = -1 + b_ + indMin;
        end
        
        period_ = positionsMin(indMin);
        period_ = period_ + 1;
        
        %Parabola values and coefficients being used for extracting the
        %period value
        if dist_(indMin) < dist_(indMin+1)
            if period_ < length(cumM_)
                if dist_(indMin) < dist_(indMin-1)
                    if period_ > 2 
                        
                        firstX = dist_(-1. + period_);
                        secX = dist_(period_);
                        thirdX = dist_(1. + period_);
                        a_= .5 * (firstX + thirdX - 2. * secX);
                        b_= .5 * (thirdX - firstX);
                        sh_= -b_ ./ a_ ./ 2; 
                        output = secX - b_ .^ 2 ./ a_ ./ 4; 
                        period_ = period_+sh_;
                        period_ = period_ - 1.;
                    end
                end
            end
        end
        
        yin.period_(k) = period_;
        
        k = k + 1;
    end

    
    fileMv(8:end)
    %Fundamental frequencies getting printed below
    fundFreqs = fs ./ yin.period_
    
    
    
    figure(u)
    l = length(fundFreqs);
    plot(timeData / l:timeData / l:timeData, fundFreqs,'-');
    hold on;
    plot(timeData / l:timeData / l:timeData, fundFreqs, 'xr');
    tit_ = strcat(fileMv(8:end), ' - my version (fund. freq.)');
    title(tit_);
    xlabel('Time (sec)');
    ylabel('Fund. freq (Hz)');
    fundFreqsTxt = num2str(floor(fundFreqs)');

    hold on;
    text(timeData / l:timeData / l:timeData, fundFreqs, fundFreqsTxt, 'VerticalAlignment','bottom', ...
        'HorizontalAlignment','right')
end


%All the vowels are stored in the following directory called 'vowels'
files = dir('pluginVals');

fileN = [];

%We iterate over the files in the directory through the below loop
for file = files';
    
    if strcmp(file.name, '.') == 0 && strcmp(file.name, '..') == 0
        fileN = [fileN; char(strcat(strcat('pluginVals\', char(file.name))))];
    end
    
end

allData = cellstr(fileN);

for k = 1:length(allData)
    
    figure(k + u);
    fileMv = char(allData(k));
    
    
    plug = dlmread(fileMv);

    timeData = plug(:, 1);
    fundFreqs = plug(:, 2);
    plot(timeData, fundFreqs,'-');
    hold on;
    plot(timeData, fundFreqs, 'xr');
        tit_ = strcat(fileMv(12:end), ' - Plug-in (fund. freq.)');
    title(tit_);
    
    xlabel('Time (sec)');
    ylabel('Fund. freq (Hz)');
    fundFreqsTxt = num2str(floor(fundFreqs));
    text(timeData, fundFreqs, fundFreqsTxt, 'VerticalAlignment','bottom', ...
                             'HorizontalAlignment','right')
end
