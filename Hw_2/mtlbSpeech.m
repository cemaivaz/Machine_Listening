    
    % Read the data back into MATLAB, and listen to audio.
    
    [y, Fs, nbits, readinfo] = wavread('u.wav');
    
    sound(y, Fs);
    
    %Below, spectral density is implemented
    
    len = 100;
    
    noverlap = 90;
    
    NFFT = 128;
    
    [y_,freq,time,p] = spectrogram(y,len,noverlap,NFFT,Fs);
    
    x = y;
    
    x1 = x.*hamming(length(x));
    
    preemph = [1 0.63];
    x1 = filter(1,preemph,x1);
    
    
    A = lpc(x1,8);
    rts = roots(A);
    
    
    rts = rts(imag(rts)>=0);
    angz = atan2(imag(rts),real(rts));
    
    [frqs,indices] = sort(angz.*(Fs/(2*pi)));
    bw = -1/2*(Fs/(2*pi))*log(abs(rts(indices)));
    
    formants = [];
    nn = 1;
    for kk = 1:length(frqs)
        if (frqs(kk) > 90 && bw(kk) <400)
            formants(nn) = frqs(kk);
            nn = nn+1;
        end
    end
    formants