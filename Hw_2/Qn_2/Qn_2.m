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
    thre_=0.11;
    maxprd=256;

    wsize=maxprd;
    skip=wsize;

    dataSize=length(x);
    x=x(:);
    
    
    frameNo=floor((dataSize-maxprd-wsize)/skip);
    pwr=ones(1,frameNo);
    prd=ones(1,frameNo);
    ap=ones(1,frameNo);
    
    % shifted data
    x=convmtx(x,maxprd+1);
    x=x(maxprd:end-maxprd,:);
    
    
    k = 1;
    while k <= frameNo
        
        st_=(k-1)*skip; % offset of frame
        xx=x(st_+1:st_+wsize,:);
        d=mean( (xx - repmat(xx(:,1),1,maxprd+1)).^2 )/2;     % squared difference function
        dd= d(2:end) ./ (cumsum(d(2:end)) ./ (1:(maxprd)));   % cumulative mean - normalized
        
        % parabolic interpolation of all triplets to refine local minima
        min_pos=1:numel(dd);    % nominal position of each sample
        x1=dd(1:end-2);
        x2=dd(2:end-1);
        x3=dd(3:end);
        a=(x1+x3-2*x2)/2;
        b=(x3-x1)/2;
        shift=-b./(2*a);        % offset of interpolated minimum re current sample
        val=x2-b.^2./(4*a);     % value of interpolated minimum
        
        % replace all local minima by their interpolated value,
        idx= 1 + find(x2<x1 & x2<x3);
        dd(idx)=val(idx-1);
        min_pos(idx)=min_pos(idx-1)+shift(idx-1);
        
        % find index of first min below thre_old
        a=dd<thre_;
        if isempty(find(a))
            [~,prd0]=min(dd); % none below thre_old, take global min instead
        else
            b=min(find(a)); % left edge
            c=min(b*2,numel(a));
            [~,prd0]=min(dd(b:(c-1)));
            prd0=b+prd0-1;
        end
        
        prd=min_pos(prd0)+1;
        
        if prd>2 & prd<numel(dd) & d(prd0)<d(prd0-1) & d(prd0)<d(prd0+1)
            
            % refine by parabolic interpolation of raw difference function
            x1=d(prd-1);
            x2=d(prd);
            x3=d(prd+1);
            a=(x1+x3-2*x2)/2;
            b=(x3-x1)/2;
            shift=-b./(2*a);        % offset of interpolated minimum re current sample
            val=x2-b.^2./(4*a);     % value of interpolated minimum
            prd=prd+shift-1;
        end
        
        % aperiodicity
        frac=prd-floor(prd);
        if frac==0
            yy=xx(:,prd);
        else
            yy=(1-frac)*xx(:,floor(prd+1))+frac*xx(:,floor(prd+1)+1); % linear interpolation
        end
        pwr=(mean(xx(:,1).^2) + mean(yy.^2))/2; % average power over fixed and shifted windows
        res=mean(((xx(:,1) - yy)).^2) / 2;
        ap=res/pwr;
        
        
        r.prd(k)=prd;
        r.ap(k)=ap;
        
        k = k + 1;
    end
    
    if nargout==0;
        subplot 211; plot(r.prd); title('period'); xlabel('frame'); ylabel('samples');
        subplot 212; plot(r.ap); title('periodicity'); xlabel('frame');
        r=[];
    end
    
    
    
    fileMv(8:end)
    fs ./ r.prd
    
end
