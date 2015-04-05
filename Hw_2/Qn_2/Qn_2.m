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
    [x fs] = wavread(fileMv);
    
    % YIN2: a simple implementation of the yin period-estimation algorithm
    %
    %  yin2(x) : plot the period, power, and aperiodicity as a function of time
    %
    %  r=yin2(x,p): use parameters in p, return result in r:
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
    
    dataSize = length(x);
    
    thre_ = 0.11;
    range_ = 256;
    
    windLen=range_;
    skip=windLen;
    
    x=x(:);
    
    
    frameNo = floor((dataSize - range_ - windLen) / skip);
    
    
    aperiod = ones(1, frameNo);
    period_ = ones(1, frameNo);
    p_ = ones(1, frameNo); %pwr
    
    % shifted data
    x = convmtx(x,range_+1);
    x = x(range_:end-range_,:);
    
    k = 1;
    while k <= frameNo
        
        st_ = k * skip - skip; % offset of frame
        xSub = x(st_ + 1:st_ + windLen,:);
        dist_ = .5 * mean(power(xSub - repmat(xSub(:,1),1,range_+1), 2));     % squared difference function
        cumM_ = dist_(2:end) ./ (cumsum(dist_(2:end)) ./ (1:(range_)));   % cumulative mean - normalized
        
        % parabolic interpolation of all triplets to refine local minima
        positionsMin=1:length(cumM_);    % nominal position of each sample
        
        firstX = cumM_(1:length(cumM_)-2);
        secX = cumM_(2:length(cumM_)-1);
        thirdX = cumM_(3:length(cumM_));
        a_ = .5 * (firstX + thirdX - 2*secX);
        b_ = .5 * (thirdX - firstX);
        sh_ = -b_./(2*a_);        % offset of interpolated minimum re current sample
        output = secX-b_.^2./(4*a_);     % value of interpolated minimum
        
        % replace all local minima by their interpolated value,
        ind_ = 1 + find(secX<firstX & secX<thirdX);
        cumM_(ind_) = output(ind_ - 1);
        positionsMin(ind_) = positionsMin(ind_ - 1) + sh_(ind_ - 1);
        
        % find index of first min below thre_old
        a_ = cumM_ < thre_;
        
        if isempty(find(a_))
            [~, indMin] = min(cumM_); % none below thre_old, take global min instead
        else
            b_ = min(find(a_)); % left edge
            c_ = min(b_*2, length(a_));
            [~, indMin] = min(cumM_(b_:(c_-1)));
            indMin=b_+indMin-1;
        end
        
        period_=positionsMin(indMin)+1;
        
        if period_>2
            if period_ < numel(cumM_)
                if dist_(indMin)<dist_(indMin-1)
                    if dist_(indMin)<dist_(indMin+1)
                        % refine by parabolic interpolation of raw difference function
                        firstX=dist_(period_-1);
                        secX=dist_(period_);
                        thirdX=dist_(period_+1);
                        a_=(firstX+thirdX-2*secX)/2;
                        b_=(thirdX-firstX)/2;
                        sh_=-b_./(2*a_);        % offset of interpolated minimum re current sample
                        output=secX-b_.^2./(4*a_);     % value of interpolated minimum
                        period_=period_+sh_-1;
                    end
                end
            end
        end
        
        % aperiodicity
        frac=period_-floor(period_);
        if frac == 0
            yy = xSub(:,period_);
        else
            yy = (1-frac)*xSub(:,floor(period_+1)) + frac * xSub(:,floor(period_ + 1)+1); % linear interpolation
        end
        p_=(mean(power(xSub(:,1), 2)) + mean(yy.^2))/2; % average power over fixed and shifted windows
        result = mean(power(((xSub(:,1) - yy)), 2)) / 2;
        aperiod = result/p_;
        
        
        yin.period_(k)=period_;
        yin.aperiod(k)=aperiod;
        
        k = k + 1;
    end
    
    
    
    fileMv(8:end)
    fs ./ yin.period_
    
end
