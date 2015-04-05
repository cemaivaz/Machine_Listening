%Cem Rýfký Aydýn 2013800054
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
    
    % YIN2: a simple implementation of the yin period-estimation algorithm
    %
    %  yin2(data) : plot the period, power, and aperiodicity as a function of time
    %
    %  r=yin2(data,p): use parameters in p, return result in r:
    %
    %    yin.prd: period
    %    yin.ap: aperiodicity measure
    %
    %    maxprd: samples, maximum of search range [default 100]
    %    minprd: samples, minimum of search range [default 2]
    %    wsize: samples, window size [default maxprd]
    %    skip: samples, frame period [default wsize]
    %    thre_: thre_old for period minimum [default: 0.1]
    %    smooth: samples, size of low-pass smoothing window [default: minprd/2]
    
    
    % defaults
    
    dataSize = length(data);
    
    thre_ = 0.11;
    range_ = 256;
    
    windLen = range_;
    skip = windLen;
    
    
    frameNo = floor((dataSize - range_ - windLen) / skip);
    
    
    aperiod = ones(1, frameNo);
    period_ = ones(1, frameNo);
    p_ = ones(1, frameNo); %pwr
    
    
    data = convmtx(data, range_ + 1.); % Data gets shifted
    
    k = 1;
    data = data(range_:end-range_,:);

    while k <= frameNo
        
        st_ = k * skip - skip; % offset of frame
        xSub = data(st_ + 1:st_ + windLen,:);
        repAmount = range_ + 1;
        dist_ = .5 * mean(power(xSub - repmat(xSub(:,1),1,repAmount), 2));     % squared difference function
        cumM_ = dist_(2:length(dist_)) ./ (cumsum(dist_(2:length(dist_))) ./ (1:(repAmount - 1)));   % cumulative mean - normalized
        
        % parabolic interpolation of all triplets to refine local minima
        positionsMin=1:length(cumM_);    % nominal position of each sample
        
        firstX = cumM_(1:length(cumM_)-2);
        secX = cumM_(2:length(cumM_)-1);
        thirdX = cumM_(3:length(cumM_));
        a_ = .5 * (firstX + thirdX - 2. * secX);
        b_ = .5 * (thirdX - firstX);
        sh_ = -b_ ./ a_ * .5;        % offset of interpolated minimum re current sample
        output = secX - b_ .^2 ./ a_ .* 0.25;     % value of interpolated minimum
        
        % replace all local minima by their interpolated value,
        ind_ = find(secX<firstX & secX<thirdX) + 1.;
        cumM_(ind_) = output(-1. + ind_);
        positionsMin(ind_) = positionsMin(-1. + ind_) + sh_(-1. + ind_);
        
        % find index of first min below thre_old
        
        a_ = zeros(length(cumM_));
        for o_ = 1:length(cumM_)
            a_(o_) = (cumM_(o_) < thre_);
        end
        
        if length(find(a_)) == 0
            [minVal, indMin] = min(cumM_); % none below thre_old, take global min instead
        else
            bTmp = find(a_);
            bTmp = bTmp(1);
            b_ = min(bTmp); % left edge
            c_ = min(b_*2., length(a_));
            [minVal, indMin] = min(cumM_(b_:(c_-1)));
            indMin = -1 + b_ + indMin;
        end
        
        period_ = positionsMin(indMin);
        period_ = period_ + 1;
        
        if dist_(indMin) < dist_(indMin+1)
            if period_ < length(cumM_)
                if dist_(indMin) < dist_(indMin-1)
                    if period_ > 2 
                        % refine by parabolic interpolation of raw difference function
                        firstX = dist_(-1. + period_);
                        secX = dist_(period_);
                        thirdX = dist_(1. + period_);
                        a_= .5 * (firstX + thirdX - 2. * secX);
                        b_= .5 * (thirdX - firstX);
                        sh_= -b_ ./ a_ ./ 2;        % offset of interpolated minimum re current sample
                        output = secX - b_ .^ 2 ./ a_ ./ 4;     % value of interpolated minimum
                        period_ = period_+sh_;
                        period_ = period_ - 1.;
                    end
                end
            end
        end
        
        % aperiodicity
        fr_ = period_;
        fr_ = fr_ - floor(period_);
        if fr_ == .0
            xRes = xSub(:,period_);
        else
            xRes = xSub(:, floor(period_ + 1)) .* (1 - fr_);
            xRes = xRes + xSub(:, floor(period_ + 1) + 1) .* fr_;% linear interpolation
        end
        p_ = .5 * (mean(power(xSub(:, 1), 2)) + mean(power(xRes, 2))); % average power over fixed and shifted windows
        result = .5 * mean(power(((xSub(:, 1) - xRes)), 2));
       
        
        yin.period_(k) = period_;
        
        k = k + 1;
    end

    fileMv(8:end)
    fs ./ yin.period_
    
end
