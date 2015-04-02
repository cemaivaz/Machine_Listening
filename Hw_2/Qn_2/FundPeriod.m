%testing a real sound , matlab code
%x=[1 3  7  2 1  9  3 1 8 ],
[x,fs,nbits]=wavread('vowels\e.wav');
sound(x,fs)%fs=44100Hz,
fs %sampling freuqncy 
start=1; %pitch a fram around t=10000
length=512;
x=x(start:start+length);
auto_corr_x=xcorr(x); %auto-correlation
figure(1), clf
subplot(2,1,1),plot(x)
title(' one frame of the sound A5-flute=880Hz')
grid on, grid(gca,'minor'), hold on
subplot(2,1,2),plot(auto_corr_x)
title('cross correlation result')
grid on, grid(gca,'minor')
 [pks,locs] = findpeaks(auto_corr_x)
[mm,peak1_ind]=max(pks)
%'peak value1 at location'
pks(peak1_ind) %peak
locs (peak1_ind) %location
 %'peak value2 at location'
pks(peak1_ind+1)%peask next to the top peak
locs (peak1_ind+1) %location
period=locs(peak1_ind+1)-locs(peak1_ind) 
 pitch_Hz=fs/period %display pitch in Hz