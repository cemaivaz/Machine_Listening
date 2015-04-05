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
    %    r.prd: period
    %    r.ap: aperiodicity measure
    %
    %    maxprd: samples, maximum of search range [default 100]
    %    minprd: samples, minimum of search range [default 2]
    %    wsize: samples, window size [default maxprd]
    %    skip: samples, frame period [default wsize]
    %    thre_: thre_old for period minimum [default: 0.1]
    %    smooth: samples, size of low-pass smoothing window [default: minprd/2]
    
    
    % defaults
    
    dataSize = length(x);
    
    thre_=0.11;
    range_=256;
    
    wsize=range_;
    skip=wsize;
    
    x=x(:);
    
    
    frameNo=floor((dataSize-range_-wsize) / skip);
    pwr = ones(1, frameNo);
    prd = ones(1, frameNo);
    ap = ones(1, frameNo);
    
    % shifted data
    x=convmtx(x,range_+1);
    x=x(range_:end-range_,:);
    
    
    k = 1;
    while k <= frameNo
        
        st_ = (k-1) * skip; % offset of frame
        xSub = x(st_+1:st_+wsize,:);
        dist_ = mean(power(xSub - repmat(xSub(:,1),1,range_+1), 2))/2;     % squared difference function
        cumM_ = dist_(2:end) ./ (cumsum(dist_(2:end)) ./ (1:(range_)));   % cumulative mean - normalized
        
        % parabolic interpolation of all triplets to refine local minima
        min_pos=1:length(cumM_);    % nominal position of each sample
        
        x1 = cumM_(1:end-2);
        x2 = cumM_(2:end-1);
        x3 = cumM_(3:end);
        a = .5 * (x1+x3-2*x2);
        b = .5 * (x3-x1);
        shift = -b./(2*a);        % offset of interpolated minimum re current sample
        val = x2-b.^2./(4*a);     % value of interpolated minimum
        
        % replace all local minima by their interpolated value,
        idx= 1 + find(x2<x1 & x2<x3);
        cumM_(idx)=val(idx-1);
        min_pos(idx)=min_pos(idx-1)+shift(idx-1);
        
        % find index of first min below thre_old
        a=cumM_<thre_;
        if isempty(find(a))
            [~,prd0]=min(cumM_); % none below thre_old, take global min instead
        else
            b=min(find(a)); % left edge
            c=min(b*2,numel(a));
            [~,prd0]=min(cumM_(b:(c-1)));
            prd0=b+prd0-1;
        end
        
        prd=min_pos(prd0)+1;
        
        if prd>2
            if prd < numel(cumM_)
                if dist_(prd0)<dist_(prd0-1)
                    if dist_(prd0)<dist_(prd0+1)
                        % refine by parabolic interpolation of raw difference function
                        x1=dist_(prd-1);
                        x2=dist_(prd);
                        x3=dist_(prd+1);
                        a=(x1+x3-2*x2)/2;
                        b=(x3-x1)/2;
                        shift=-b./(2*a);        % offset of interpolated minimum re current sample
                        val=x2-b.^2./(4*a);     % value of interpolated minimum
                        prd=prd+shift-1;
                    end
                end
            end
        end
        
        % aperiodicity
        frac=prd-floor(prd);
        if frac==0
            yy=xSub(:,prd);
        else
            yy=(1-frac)*xSub(:,floor(prd+1))+frac*xSub(:,floor(prd+1)+1); % linear interpolation
        end
        pwr=(mean(xSub(:,1).^2) + mean(yy.^2))/2; % average power over fixed and shifted windows
        res=mean(((xSub(:,1) - yy)).^2) / 2;
        ap=res/pwr;
        
        
        r.prd(k)=prd;
        r.ap(k)=ap;
        
        k = k + 1;
    end
    
    
    
    fileMv(8:end)
    fs ./ r.prd
    
end
