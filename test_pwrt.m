fs = 10000;
t = 0:1/fs:0.3-1/fs;

l = [0 130.81 146.83 164.81 174.61 196.00 220 246.94];
m = [0 261.63 293.66 329.63 349.23 392.00 440 493.88];
h = [0 523.25 587.33 659.25 698.46 783.99 880 987.77];
note = @(f,g) [1 1 1]*sin(2*pi*[l(g) m(g) h(f)]'.*t);

mel = [3 2 1 2 3 3 3 0 2 2 2 0 3 5 5 0 3 2 1 2 3 3 3 3 2 2 3 2 1]+1;
acc = [3 0 5 0 3 0 3 3 2 0 2 2 3 0 5 5 3 0 5 0 3 3 3 0 2 2 3 0 1]+1;

song = [];
for kj = 1:length(mel)
    song = [song note(mel(kj),acc(kj)) zeros(1,0.01*fs)];
end
song = song/(max(abs(song))+0.1);

% To hear, type sound(song,fs)

x = pspectrum(b_sum,10000,'spectrogram','TimeResolution',0.1, ...
    'OverlapPercent',0,'MinThreshold',-10);


y = bandpass(b_sum,[100 2000], 10000);

pspectrum(y,10000,'spectrogram','TimeResolution',0.1, ...
    'OverlapPercent',0,'MinThreshold',-1)

y = bandpass(b_sum,[100 2000], 10000);
y =  notch(y, 10000, 1000,1);
 
 pspectrum(y,10000,'spectrogram','TimeResolution',0.1, ...
    'OverlapPercent',0,'MinThreshold',-20)



soundsc(y(1:50000))


hong = lowpass(b_sum,1000,fs);
hong = highpass(hong,100,fs);

pspectrum(hong,10000,'spectrogram','TimeResolution',0.1, ...
    'OverlapPercent',0,'MinThreshold',-10)



      

hong_zero_pad = hong;
hong_zero_pad(hong_zero_pad<0) = 0

pspectrum(hong_zero_pad,10000,'spectrogram','TimeResolution',0.1, ...
    'OverlapPercent',0,'MinThreshold',-10)


[P,F,T] = pspectrum(hong,10000,'spectrogram','TimeResolution',0.1, ...
    'OverlapPercent',0,'MinThreshold',-50);

soundsc(hong_zero_pad(1:10000)*1000)



ifft()



pspectrum(b,10000,'spectrogram','TimeResolution',0.1,'OverlapPercent',0,'MinThreshold',-60);
soundsc(b(1:10000))
b = lowpass(b,100,fs);
soundsc(b)


subplot(5,1,1)
pspectrum(b_sum,10000,'spectrogram','TimeResolution',0.1,'OverlapPercent',0,'MinThreshold',-20);
set(gca,'fontsize',16)
title('IPS')



soundsc(b_sum(10000:50000))


y = b_sum;
y = bandpass(b_sum,[300 2000], 10000);
soundsc(y(20000:60000))

