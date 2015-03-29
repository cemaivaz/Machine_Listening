
function [coef]=comp_dgt_fb(f,g,a,M)
%COMP_DGT_FB  Filter bank DGT
%   Usage:  c=comp_dgt_fb(f,g,a,M);
%  
%   This is a computational routine. Do not call it directly.
%
%   See help on DGT.

%   AUTHOR : Peter L. S?¸ndergaard.

% Calculate the parameters that was not specified.

%{
a = a

M = M
%}

L=size(f,1);

N=L/a;
gl=length(g);
W=size(f,2);      % Number of columns to apply the transform to.;
glh=floor(gl/2);  % gl-half



% Conjugate the window here.
g=conj(fftshift(g));

coef=zeros(M,N,W,'double');

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