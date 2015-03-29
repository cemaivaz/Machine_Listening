function [pos,V0] = onsetdet(f,win_length,thre,range,multi,shift,showplot)
%ONSETDET  Onset detection wrapper
%   Usage: [pos,V0] = onsetdet(f,win_length,thre,range,multi,shift,showplot)
%          [pos,V0] = onsetdet(f,win_length,thre,range,multi,shift)
%          [pos,V0] = onsetdet(f,win_length,thre,range,multi)
%          [pos,V0] = onsetdet(f,win_length,thre,range)
%          [pos,V0] = onsetdet(f,win_length,thre)
%          [pos,V0] = onsetdet(f,win_length)
%          pos = onsetdet(...)
%
%   Input parameters:
%         f         : Signal to be analyzed (single channel only)
%         win_length: Window length for the STFT analysis (in samples)
%         thre      : Peak-picking threshold
%         range     : Area of interest for the choice of local maxima
%         multi     : Area of interest multiplication factor for the
%                     peak-picking
%         shift     : Readjustment of the peaks
%                     (in shift~cdot~win_length/16)
%         showplot  : Plot the results (0/1)
%   Output parameters:
%         pos       : Onset sequence
%         V0        : Regular discrete Gabor transform of f
%
%   This routine produces a sequence of onsets using a straightforward
%   realization of a spectral flux based onset detection process as
%   described, e.g. by Dixon (see reference).
%
%   The spectral flux onset detection function is computed with a 16 times
%   redundant Gabor transform using a Hann window, implemented in
%   SPECFLUX.
%
%   Local maxima of the onset detection function are chosen as onset if
%   larger than the local mean by at least the threshold parameter thre.
%   This choice is performed by PEAKPICK, a simple peakpicking algorithm.
%
%   A time slice is considered a local maximum if its spectral flux value
%   is larger than those of the surrounding slices on an area of +-range.
%   The local mean is computed as the mean value of the spectral flux
%   function on an area corresponding to -multi*range to +range of the
%   current position.
%
%   See also:  onsetnsgt, invonsetnsgt, specflux, peakpick
%
%   References:
%     S. Dixon. Onset detection revisited. In Proceedings of the 9th
%     International Conference on Digital Audio Effects, volume 120, pages
%     133-137, 2006.
%
%     P. Balazs, M. Dörfler, F. Jaillet, N. Holighaus, and G. A. Velasco.
%     Theory, implementation and applications of nonstationary Gabor Frames.
%     J. Comput. Appl. Math., 236(6):1481-1496, 2011.
%
%
%   Url: http://nsg.sourceforge.net/doc/helpers/onsetdet.php

% Copyright (C) 2013 Nicki Holighaus.
% This file is part of NSGToolbox version 0.1.0
%
% This work is licensed under the Creative Commons
% Attribution-NonCommercial-ShareAlike 3.0 Unported
% License. To view a copy of this license, visit
% http://creativecommons.org/licenses/by-nc-sa/3.0/
% or send a letter to
% Creative Commons, 444 Castro Street, Suite 900,
% Mountain View, California, 94041, USA.

% Author: Nicki Holighaus
% Date: 26.04.13

% Check input arguments and set default arguments if necessary

if nargin < 7
    showplot = 0;
    if nargin < 6
        shift = 0;
        if nargin < 5
            multi = 3;
            if nargin < 4
                range = 3;
                if nargin < 3
                    error('Not enough input arguments');
                end
            end
        end
    end
end

[Ls,col] = size(f);

if min(Ls,col) > 1
    error('Right now, this routine supports only single channel signals');
end

if ( col ~= 1 && Ls == 1 )
    f = f.';
    Ls = col;
end

clear col;

% Compute the spectral flux

tgap = win_length/16;
[ODF,V0]=specflux(f,win_length,tgap);

% Select the significant paeks in the spectral
% flux

pos = peakpick(ODF,thre,range,multi);

% Shift the onset positions according by a fixed amount to be more precise
% (experimental, but improves the results on simple signals)

pos = onsets(pos,tgap,win_length,shift);

% Due to periodization and shifts, some onsets might appear after the
% end of the signal. Those are omitted

X = length(pos);
if ( X > 0 )
    while ( pos(X) >= Ls)
        X = X-1;
    end
end

% The first sample is always considered an onset

pos = [1,pos(1:X)].';

if showplot ~= 0    % Plot the results
    
    pos2 = floor(1+(pos-1)/tgap);
    g=max(ODF)*ones(length(ODF),1);
    g(pos2) = min(ODF)*ones(length(pos2),1);
    
    subplot(1,2,1);
    imagesc(20*log10(abs(V0(size(V0,1)*4/8+1:size(V0,1),:))+eps));
    vline(pos2,'k');
    subplot(1,2,2);
    plot(ODF,'k'); hold on; plot(g,'r+');
    axis tight; hold off; shg
end

end

% onsets should shift the detected onsets uniformly
% to another position. In many cases, this might improve
% results significantly.
%
% The real onsets do not always exactly coincide with
% the chosen peaks, but for simple signals they should
% be around (peak + some constant times [+-tgap]).
%
% Right now, this is just an experimental routine
% and still a work in progress. Work has to be done.

function pos = onsets(peaks,tgap,win_length,shift)

pos = 1+(peaks-1).*tgap + floor(shift*win_length/16);

end


function [SF,V0] = specflux(f,win_length,tgap)
%SPECFLUX  Spectral flux onset detection function
%   Usage: [SF,V0] = specflux(f,win_length,tgap)
%          SF = specflux(f,win_length,tgap)
%
%   Input parameters: 
%         f         : Input signal
%         win_length: Desired window length for the STFT
%         tgap      : Time step for the STFT
%   Output parameters:
%         SF        : Spectral flux of f*
%         V0        : STFT coefficients of f*
% 
%   This is a helper function for ONSETDET and not meant to
%   be used individually.
%
%   Computes the spectral flux onset-detection function
%   of f with a Hann window of length win_length. 
%   The STFT is taken with time shift parameter tgap*
%   and win_length frequency channels.
%
%   Externals: COMP_DGT_FB (LTFAT routine, included in NSGToolbox V0.1.0 
%              and higher)
%
%   See also:  onsetdet
%
%   References:
%     S. Dixon. Onset detection revisited. In Proceedings of the 9th
%     International Conference on Digital Audio Effects, volume 120, pages
%     133-137, 2006.
%     
%
%   Url: http://nsg.sourceforge.net/doc/helpers/specflux.php

% Copyright (C) 2013 Nicki Holighaus.
% This file is part of NSGToolbox version 0.1.0
% 
% This work is licensed under the Creative Commons 
% Attribution-NonCommercial-ShareAlike 3.0 Unported 
% License. To view a copy of this license, visit 
% http://creativecommons.org/licenses/by-nc-sa/3.0/ 
% or send a letter to 
% Creative Commons, 444 Castro Street, Suite 900, 
% Mountain View, California, 94041, USA.

% Author: Nicki Holighaus
% Date: 26.04.13

% Check input arguments

if nargin < 3
    error('Not enough input arguments');
end

% Compute the Gabor transform (sampled STFT) of f

win=winfuns('hann',win_length);

if size(f,2) > 1
    f = f.';
end

V0 = comp_dgt_fb(f,win,tgap,win_length);

% Compute the spectral flux 

VV = abs(V0);
VV = max(VV-circshift(VV,[0,1]),0);

SF = sum(VV);

% Normalize

SF = SF-mean(SF);
SF = SF./std(SF);

end

function peaks = peakpick(SF,thre,range,multi)
%PEAKPICK  Peakpicking routine
%   Usage:  peaks = peakpick(SF,thre,range,multi)
%
%   Input parameters: 
%         SF        : Onset detection function
%         thre      : Threshold value
%         range     : Relevance area for local maximum
%         multi     : Asymmetric extension factor for relevance area
%   Output parameters:
%         peaks     : Significant maxima of SF*
%
%   This is a helper function for ONSETDET and not meant to be used 
%   individually.
%
%   For an onset detection function SF, the routine picks only those 
%   local maxima that are larger than the local mean over an area of the 
%   form
%
%           <---multi*range---X---range--->
%
%   by more than the threshold given by thre.
%
%   See also:  onsetdet
%
%   References:
%     S. Dixon. Onset detection revisited. In Proceedings of the 9th
%     International Conference on Digital Audio Effects, volume 120, pages
%     133-137, 2006.
%     
%
%   Url: http://nsg.sourceforge.net/doc/helpers/peakpick.php

% Copyright (C) 2013 Nicki Holighaus.
% This file is part of NSGToolbox version 0.1.0
% 
% This work is licensed under the Creative Commons 
% Attribution-NonCommercial-ShareAlike 3.0 Unported 
% License. To view a copy of this license, visit 
% http://creativecommons.org/licenses/by-nc-sa/3.0/ 
% or send a letter to 
% Creative Commons, 444 Castro Street, Suite 900, 
% Mountain View, California, 94041, USA.

% Author: Nicki Holighaus
% Date: 26.04.13

% Check input arguments

if nargin < 4
    error('Not enough input arguments');
end

% Compute local maxima

maxima = locmax(SF,range);

m = length(maxima);
n = length(SF);

% Since the signal is considered periodic, some
% periodic border values are added

SF = [SF(n-(multi+1)*range+1:1:n),SF,SF(1:range)];

kk=0;
peaks = [];

% The loop selects the significant local maxima
% from the output of locmax

for ii = 1:m
    pos = maxima(ii);
    th_loc_mean = sum(SF(pos:pos+(multi+2)*range))/((multi+2)*range+1)+thre;
    if ( SF(pos+(multi+1)*range) > th_loc_mean )
        kk = kk + 1;
        peaks(kk) = maxima(ii);
    end
end   

end

% Computes 'some kind of' local maxima.
% A point of SF is considered a local maximum
% in the sense of this routine, if it is larger
% than all surrounding points in a neigborhood
% with a radius of range-1 of this point.

function maxima = locmax(SF,range)

if ( size(SF,1) ~= 1 )
    SF=SF.';
end

n = length(SF);

SF = [SF(n-range+1:1:n),SF,SF(1:range)];

kk = 1;
maxima = [];

for ii = 1:n
    xx = SF(ii:ii+2*range);
    yy = ones(1,2*range+1)*SF(ii+range);
    if ( sum( yy >= xx ) == 2*range+1 )
        maxima(kk) = ii; 
        kk = kk+1;
    end   
end

end

function [coef]=comp_dgt_fb(f,g,a,M)
%COMP_DGT_FB  Filter bank DGT
%   Usage:  c=comp_dgt_fb(f,g,a,M);
%  
%   This is a computational routine. Do not call it directly.
%
%   See help on DGT.

%   AUTHOR : Peter L. S?¸ndergaard.

% Calculate the parameters that was not specified.
L=size(f,1);
N=L/a;
gl=length(g);
W=size(f,2);      % Number of columns to apply the transform to.
glh=floor(gl/2);  % gl-half


% Conjugate the window here.
g=conj(fftshift(g));

coef=zeros(M,N,W,assert_classname(f,g));

% ----- Handle the first boundary using periodic boundary conditions. ---
for n=0:ceil(glh/a)-1

    % Periodic boundary condition.
    fpart=[f(L-(glh-n*a)+1:L,:);...
           f(1:gl-(glh-n*a),:)];
    
    fg=bsxfun(@times,fpart,g);
    
    % Do the sum (decimation in frequency, Poisson summation)
    coef(:,n+1,:)=sum(reshape(fg,M,gl/M,W),2);
      
end;

% ----- Handle the middle case. ---------------------
for n=ceil(glh/a):floor((L-ceil(gl/2))/a)
  
  fg=bsxfun(@times,f(n*a-glh+1:n*a-glh+gl,:),g);
  
  % Do the sum (decimation in frequency, Poisson summation)
  coef(:,n+1,:)=sum(reshape(fg,M,gl/M,W),2);
end;

% ----- Handle the last boundary using periodic boundary conditions. ---
for n=floor((L-ceil(gl/2))/a)+1:N-1

    % Periodic boundary condition.
    fpart=[f((n*a-glh)+1:L,:);... %   L-n*a+glh elements
           f(1:n*a-glh+gl-L,:)];  %  gl-L+n*a-glh elements      
    
    fg=bsxfun(@times,fpart,g);
    
    % Do the sum (decimation in frequency, Poisson summation)
    coef(:,n+1,:)=sum(reshape(fg,M,gl/M,W),2);      
end;

% --- Shift back again to make it a frequency-invariant system. ---
for n=0:N-1
  coef(:,n+1,:)=circshift(coef(:,n+1,:),n*a-glh);
end;


coef=fft(coef);

end



% Simple code using a lot of circshifts.
% Move f initially so it lines up with the initial fftshift of the
% window
%f=circshift(f,glh);
%for n=0:N-1
  % Do the inner product.
  %fg=circshift(f,-n*a)(1:gl,:).*gw;
  
  % Periodize it.
  %fpp=zeros(M,W);
  %for ii=0:gl/M-1
    %  fpp=fpp+fg(ii*M+1:(ii+1)*M,:);
    %end;
%  fpp=sum(reshape(fg,M,gl/M,W),2);
  
  % Shift back again.
%  coef(:,n+1,:)=circshift(fpp,n*a-glh); %),M,1,W);
  
%end;



function g = winfuns(name,x,L)
%WINFUNS  Window function generator  
%   Usage:  g = winfuns(name,x)
%           g = winfuns(name,N,L)
%           g = winfuns(name,N)
%   
%   Input parameters: 
%         name      : String containing the window name
%         x         : Vector of sampling positions
%         N         : Window support (in samples)
%         L         : Output length (in samples)
%   Output parameters:
%         g         : Output window
%   
%   This function serves to compute a variety of standard and some more 
%   exotic window functions. Most of the functions used are detailed and 
%   discussed in classical papers (see references below), but several are
%   included for special purposes in the toolbox only.
%   
%   Given a character string name containing the name of the desired
%   window function, the function offers 2 modes of operation. If the 
%   second input parameter is a vector x of sampling values, then the
%   specified function is evaluated at the given points. If a window length
%   N and optionally a signal length L are supplied, a symmetric, 
%   whole-point centered window with a support of N samples is produced 
%   and, given L, zero-extended to length L.
%
%   The following windows are available:
%
%     'hann'         von Hann window. Forms a PU. The Hann window has a
%                    mainlobe with of 8/N, a PSL of -31.5 dB and decay rate
%                    of 18 dB/Octave.
%
%     'cos'          Cosine window. This is the square root of the Hanning
%                    window. The cosine window has a mainlobe width of 6/N,
%                    a  PSL of -22.3 dB and decay rate of 12 dB/Octave.
%                  
%     'rec'          Rectangular window. The rectangular window has a
%                    mainlobe width of 4/N, a  PSL of -13.3 dB and decay
%                    rate of 6 dB/Octave. Forms a PU. Alias: 'square'
%
%     'tri'          Triangular window. 
%
%     'hamming'      Hamming window. Forms a PU that sums to 1.08 instead
%                    of 1.0 as usual. The Hamming window has a
%                    mainlobe width of 8/N, a  PSL of -42.7 dB and decay
%                    rate of 6 dB/Octave.
%
%     'blackman'     Blackman window. The Blackman window has a
%                    mainlobe width of 12/N, a PSL of -58.1 dB and decay
%                    rate of 18 dB/Octave.
%
%     'blackharr'    Blackman-Harris window. The Blackman-Harris window has 
%                    a mainlobe width of 16/N, a PSL of -92.04 dB and decay
%                    rate of 6 dB/Octave.
%
%     'modblackharr'  Modified Blackman-Harris window. This slightly 
%                     modified version of the Blackman-Harris window has 
%                     a mainlobe width of 16/N, a PSL of -90.24 dB and decay
%                     rate of 18 dB/Octave.
%
%     'nuttall'      Nuttall window. The Nuttall window has a mainlobe 
%                    width of 16/N, a PSL of -93.32 dB and decay rate of 
%                    18 dB/Octave.
%
%     'nuttall10'    2-term Nuttall window with 1 continuous derivative. 
%                    Alias: 'hann'.
%
%     'nuttall01'    2-term Nuttall window with 0 continuous derivatives. 
%                    Alias: 'hamming'.
%
%     'nuttall20'    3-term Nuttall window with 3 continuous derivatives. 
%                    The window has a mainlobe width of 12/N, a PSL of 
%                    -46.74 dB and decay rate of 30 dB/Octave.
%
%     'nuttall11'    3-term Nuttall window with 1 continuous derivative. 
%                    The window has a mainlobe width of 12/N, a PSL of 
%                    -64.19 dB and decay rate of 18 dB/Octave.
%
%     'nuttall02'    3-term Nuttall window with 0 continuous derivatives. 
%                    The window has a mainlobe width of 12/N, a PSL of 
%                    -71.48 dB and decay rate of 6 dB/Octave.
%
%     'nuttall30'    4-term Nuttall window with 5 continuous derivatives. 
%                    The window has a mainlobe width of 16/N, a PSL of 
%                    -60.95 dB and decay rate of 42 dB/Octave.
%
%     'nuttall21'    4-term Nuttall window with 3 continuous derivatives. 
%                    The window has a mainlobe width of 16/N, a PSL of 
%                    -82.60 dB and decay rate of 30 dB/Octave.
%
%     'nuttall12'    4-term Nuttall window with 1 continuous derivatives. 
%                    Alias: 'nuttall'.
%
%     'nuttall03'    4-term Nuttall window with 0 continuous derivatives. 
%                    The window has a mainlobe width of 16/N, a PSL of 
%                    -98.17 dB and decay rate of 6 dB/Octave.
%
%     'gauss'        Truncated, stretched Gaussian: exp(-18*x^2) restricted
%                    to the interval ]-.5,.5[.
%
%     'wp2inp'       Warped Wavelet uncertainty equalizer (see WP 2 of the
%                    EU funded project UnlocX). This function is included 
%                    as a test function for the Wavelet transform 
%                    implementation and serves no other purpose in this 
%                    toolbox.
%
%   See also:  nsgcqwin, nsgwvltwin, nsgerbwin
%
%   References:
%     Wikipedia. Window function - wikipedia article.
%     http://en.wikipedia.org/wiki/Window_function.
%     
%     A. Nuttall. Some windows with very good sidelobe behavior. IEEE Trans.
%     Acoust. Speech Signal Process., 29(1):84-91, 1981.
%     
%     F. Harris. On the use of windows for harmonic analysis with the
%     discrete Fourier transform. Proceedings of the IEEE, 66(1):51 - 83,
%     January 1978.
%     
%
%   Url: http://nsg.sourceforge.net/doc/windows/winfuns.php

% Copyright (C) 2013 Nicki Holighaus.
% This file is part of NSGToolbox version 0.1.0
% 
% This work is licensed under the Creative Commons 
% Attribution-NonCommercial-ShareAlike 3.0 Unported 
% License. To view a copy of this license, visit 
% http://creativecommons.org/licenses/by-nc-sa/3.0/ 
% or send a letter to 
% Creative Commons, 444 Castro Street, Suite 900, 
% Mountain View, California, 94041, USA.

% Author: Nicki Holighaus
% Date: 25.04.13

if nargin < 2
    error('Not enough input arguments');
end

if numel(x) == 1
    N = x;
    if nargin < 3
        L = N;
    end
    if L<N
        error('Output length L must be larger than or equal to N');
    end
    if mod(N,2) == 0 % For even N the sampling interval is [-.5,.5-1/N]
        x = [0:1/N:.5-1/N,-N*ones(1,L-N),-.5:1/N:-1/N]';
    else % For odd N the sampling interval is [-.5+1/(2N),.5-1/(2N)]
        x = [0:1/N:.5-.5/N,-N*ones(1,L-N),-.5+.5/N:1/N:-1/N]';
    end
end

if size(x,2) > 1
    x = x.';
end

switch name    
    case {'Hann','hann','nuttall10','Nuttall10'}
        g = .5 + .5*cos(2*pi*x);
        
    case {'Cosine','cosine','cos','Cos','sqrthann','Sqrthann'}
        g = cos(pi*x);
        
    case {'hamming','nuttall01','Hamming','Nuttall01'}
        g = .54 + .46*cos(2*pi*x);
        
    case {'square','rec','Square','Rec'}
        g = double(abs(x) < .5);
        
    case {'tri','triangular','bartlett','Tri','Triangular','Bartlett'}
        g = 1-2*abs(x);
        
    case {'blackman','Blackman'}
        g = .42 + .5*cos(2*pi*x) + .08*cos(4*pi*x);
        
    case {'blackharr','Blackharr'}
        g = .35875 + .48829*cos(2*pi*x) + .14128*cos(4*pi*x) + ...
            .01168*cos(6*pi*x);
        
    case {'modblackharr','Modblackharr'}
        g = .35872 + .48832*cos(2*pi*x) + .14128*cos(4*pi*x) + ...
            .01168*cos(6*pi*x);
        
    case {'nuttall','nuttall12','Nuttall','Nuttall12'}
        g = .355768 + .487396*cos(2*pi*x) + .144232*cos(4*pi*x) + ...
            .012604*cos(6*pi*x);
        
    case {'nuttall20','Nuttall20'}
        g = 3/8 + 4/8*cos(2*pi*x) + 1/8*cos(4*pi*x);
        
    case {'nuttall11','Nuttall11'}
        g = .40897 + .5*cos(2*pi*x) + .09103*cos(4*pi*x);
        
    case {'nuttall02','Nuttall02'}
        g = .4243801 + .4973406*cos(2*pi*x) + .0782793*cos(4*pi*x);
        
    case {'nuttall30','Nuttall30'}
        g = 10/32 + 15/32*cos(2*pi*x) + 6/32*cos(4*pi*x) + ...
            1/32*cos(6*pi*x);
        
    case {'nuttall21','Nuttall21'}
        g = .338946 + .481973*cos(2*pi*x) + .161054*cos(4*pi*x) + ...
            .018027*cos(6*pi*x);
        
    case {'nuttall03','Nuttall03'}
        g = .3635819 + .4891775*cos(2*pi*x) + .1365995*cos(4*pi*x) + ...
            .0106411*cos(6*pi*x);
        
    case {'gauss','truncgauss','Gauss','Truncgauss'}
        g = exp(-18*x.^2);
        
    case {'wp2inp','Wp2inp'}
        g = exp(exp(-2*x)*25.*(1+2*x));
        g = g/max(g);
        
    otherwise
        error('Unknown window function: %s.',name);
end;

% Force the window to 0 outside (-.5,.5)
g = g.*(abs(x) < .5);  