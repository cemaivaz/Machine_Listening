% load mtlb
% 
% segmentlen = 100;
% noverlap = 90;
% NFFT = 128;
% 
% spectrogram(mtlb,segmentlen,noverlap,NFFT,Fs,'yaxis')
% 
% dt = 1/Fs;
% I0 = round(0.1/dt);
% Iend = round(0.25/dt);
% x = mtlb(I0:Iend);
% 
% x1 = x.*hamming(length(x));
% 
% preemph = [1 0.63];
% x1 = filter(1,preemph,x1);
% 
% 
% A = lpc(x1,8);
% rts = roots(A);
% 
% 
% rts = rts(imag(rts)>=0);
% angz = atan2(imag(rts),real(rts));
% 
% [frqs,indices] = sort(angz.*(Fs/(2*pi)));
% bw = -1/2*(Fs/(2*pi))*log(abs(rts(indices)));
% 
% nn = 1;
% for kk = 1:length(frqs)
%     if (frqs(kk) > 90 && bw(kk) <400)
%         formants(nn) = frqs(kk);
%         nn = nn+1;
%     end
% end
% formants




    load mtlb;
    dt = 1/Fs;
I0 = round(0.1/dt);
Iend = round(0.25/dt);
x = mtlb(I0:Iend);
vals = x;
    
    x1 = vals.*hamming(length(vals));
    
    %Degisik preemph'ler icin dene
    preemph = [1 0.63];
    
    x1 = filter(1,preemph,x1);
    
    
    formNo = 2;
    A = lpc(x1, formNo * 2 + 2);
    
    
    roots_ = roots(A);
    
    roots_ = roots_(imag(roots_)>=0);
    
    angz = atan2(imag(roots_),real(roots_));
    
    [freqs,index] = sort(angz.*(Fs/(2*pi)));
    
    bw = -1/2*(Fs/(2*pi))*log(abs(roots_(index)));
    
    i = 1;
    
    [c3, ind] = max(size(freqs));
    formantVals = [];
    
    for j = 1:c3
        
        if (freqs(j) > 90 && bw(j) <400)
            
            formantVals(i) = freqs(j);
            
            i = i+1;
            
        end
        
    end
    
    formantVals(end:-1:end - 1)