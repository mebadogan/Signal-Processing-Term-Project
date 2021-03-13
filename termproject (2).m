close all; clear all; clc;
%Coders: Mehmet Baris Doganoglu and Murat Sokucu
%% Signal values and spectrums
[y,Fs] = audioread("dtmf.wav");
Numsamps = length(y);
t = 0:length(y)-1;          %Preparing time data for plot
%Fourier Transform
yffty = abs(fft(y));            %Absolute magnitude to get rid of imaginary parts
yfft1 = (yffty(1:Numsamps)*(max(y)/max(yffty))); %Discard Half of Points
f = Fs*(1:Numsamps)/Numsamps;   %Calculating frequency data for plot
fc=2000;
%Plot in Time Domain
figure (1)
plot(t, y)
xlabel('Time (seconds)')
ylabel('Amplitude')
title('Signal in Time Domain')

%Plot in Frequency Domain
figure (2)
plot(f,yfft1) 
xlabel('Frequency (Hz)')
ylabel('Amplitude')
title('Signal in Frequency domain')

% Modulation
%% AM
% Time domain
yc=1*cos(2*pi*fc*t);
mmt=(10+1*y').*yc;
figure (3)
plot(t,mmt)
xlabel('Time (seconds)')
ylabel('Amplitude')
title('AM Modulated Signal in Time Domain')
%frequency domain
amfft=abs(fft(mmt));
amfft(1)=mean(amfft(2:end)); %?
amfft1=(amfft(1:Numsamps)*(max(mmt)/max(amfft)));
figure (4)
plot(f,amfft1);
ylim([0 12])
xlim([0 8000])
title('AM Modulated Signal in Frequency Domain')
%% FM
%Time Domain
B=obw(y);
kf=((2*pi*B)/max(y)); %Assuming modulation index as 1
fmt=1*cos((2*pi*fc*t)+kf*cumsum(y'));
figure (5)
plot(t,fmt);
title('FM Modulated Signal in Time Domain')
%Frequency Domain 
fmyffty = abs(fft(fmt));            %Absolute magnitude to get rid of imaginary parts
yfft2 = (fmyffty(1:Numsamps)/200);%*(max(fmt)/max(yffty)); %Discard Half of Points
f = Fs*(1:Numsamps)/Numsamps;   %Calculating frequency data for plot
figure (6)
plot(f,yfft2);
ylim([0 1]);
title('FM Modulated Signal in Frequency Domain')
%% Wireless Channel for AM
sampleRate500kHz = 500e3;    % Sample rate of 500K Hz
sampleRate20kHz  = 20e3;     % Sample rate of 20K Hz
maxDopplerShift  = 200;      % Maximum Doppler shift of diffuse components (Hz)
delayVector = (0:5:15)*1e-6; % Discrete delays of four-path channel (s)
gainVector  = [0 -3 -6 -9];  % Average path gains (dB)

%Configure a Rayleigh channel object
rayChan = comm.RayleighChannel('SampleRate',sampleRate500kHz,'PathDelays',delayVector,'AveragePathGains',gainVector,'MaximumDopplerShift', maxDopplerShift,'Seed',10,'PathGainsOutputPort', true);

release(rayChan);
rayChan.Visualization = 'Impulse and frequency responses';
rayChan.SamplesToDisplay = '100%';

numFrames =2;
for i = 1:numFrames % Display impulse and frequency responses for 2 frames
rayChan(mmt');
amrayChanOut = rayChan(mmt');
end

figure (7)
plot(t,amrayChanOut');
%% Wireless Channel for FM
sampleRate500kHz = 500e3;    % Sample rate of 500K Hz
sampleRate20kHz  = 20e3;     % Sample rate of 20K Hz
maxDopplerShift  = 200;      % Maximum Doppler shift of diffuse components (Hz)
delayVector = (0:5:15)*1e-6; % Discrete delays of four-path channel (s)
gainVector  = [0 -3 -6 -9];  % Average path gains (dB)

%Configure a Rayleigh channel object
rayChan = comm.RayleighChannel('SampleRate',sampleRate500kHz,'PathDelays',delayVector,'AveragePathGains',gainVector,'MaximumDopplerShift', maxDopplerShift,'Seed',10,'PathGainsOutputPort', true);

release(rayChan);
rayChan.Visualization = 'Impulse and frequency responses';
rayChan.SamplesToDisplay = '100%';

numFrames =2;
for i = 1:numFrames % Display impulse and frequency responses for 2 frames
rayChan(fmt');
fmrayChanOut = rayChan(fmt');
end

figure (8)
plot(t,fmrayChanOut');
%% Noise (Gaussian)
%% AM Noising
%Plot in Time Domain with noise
figure(9)
noisedam = awgn(mmt,50);
plot(t, noisedam)
xlabel('Time (seconds)')
ylabel('Amplitude')
title('AM Modulated and Noised Signal in Time Domain')

%Fourier Transform
yfftx = abs(fft(noisedam));                                %Absolute magnitude to get rid of imaginary parts
yfftx1 = yfftx(1:Numsamps)*(max(noisedam)/max(yfftx(2:end)));      %Discard Half of Points
f = Fs*(1:Numsamps)/Numsamps;                           %Calculating frequency data for plot

%Plot in Frequency Domain with noise
figure(10)
plot(f,yfftx1)
ylim([0 15])
xlim([0 8000])
xlabel('Frequency (Hz)')
ylabel('Amplitude')
title('AM Modulated and Noised Signal in Frequency domain')
%% FM Noising
figure(11)
noisedfm = awgn(fmt,100);
plot(t, noisedfm)
xlabel('Time (seconds)')
ylabel('Amplitude')
title('FM Modulated and Noised Signal in Time Domain')

%Fourier Transform
yfmftk = abs(fft(noisedfm));                                %Absolute magnitude to get rid of imaginary parts
yfmftk1 = yfmftk((1:Numsamps));      %Discard Half of Points
                       
%Plot in Frequency Domain with noise
figure (12)
plot(f,(yfmftk1/210))
ylim([0 1])
%ylim([0 15])
%xlim([0 8000])
xlabel('Frequency (Hz)')
ylabel('Amplitude')
title('FM Modulated and Noised Signal in Frequency domain')
%% AM demodulation and Filtering
dmmt=noisedam.*yc;
[b,a] = butter(1,2000/4000);   % To define boundaries of the filter                
fdmmt = filter(b,a,dmmt);       % Low filter
dfdmmt=detrend(fdmmt);
figure(13)
plot(detrend(dfdmmt)); %detrend(data) computes and subtracts the mean value from each time-domain signal in data
ylim([-1,+1]);
title('Filtered and AM Demodulated Signal in Time Domain')
damyffty = abs(fft(dfdmmt));            %Absolute magnitude to get rid of imaginary parts
damyfft = (damyffty(1:Numsamps)*(max(dfdmmt)/max(damyffty))); %Discard Half of Points
figure(14)
plot(damyfft);
title('Filtered and AM Demodulated Signal in Frequency Domain')

%% FM demodulation and filtering
dem = diff(noisedfm);                 
ffdmmt = bandpass(dem,[100 8000],Fs);
ddffdmmt=(detrend(ffdmmt)/2);
figure(15)
plot(ddffdmmt);
title('Filtered and FM Demodulated Signal in Time Domain')
dfmyffty = abs(fft(ddffdmmt));            %Absolute magnitude to get rid of imaginary parts
dfmyfft = dfmyffty(1:Numsamps-1)*(max(ddffdmmt)/max(dfmyffty)); %Discard Half of Points
figure(16)
plot(dfmyfft);
title('Filtered and FM Demodulated Signal in Frequency Domain')
dem2 = diff(noisedfm);
[up,low]=envelope(dem2)
figure(17)
plot(t(1:6399),up/2)
hold on 
plot(t(1:6399),low/2)
hold off
title('Filtered and FM Demodulated Signal in Time Domain')
fup = abs(fft(up));            %Absolute magnitude to get rid of imaginary parts
fup = fup(1:Numsamps-1); %Discard Half of Points
flow= abs(fft(low));            %Absolute magnitude to get rid of imaginary parts
flow = flow(1:Numsamps-1); %Discard Half of Points
figure(18)
plot(f(1:6399),fup/200)
ylim([0 1])
hold on 
plot(f(1:6399),flow/200)
hold off
title('Filtered and FM Demodulated Signal in Frequency Domain')