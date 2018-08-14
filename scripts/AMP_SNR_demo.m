clear all;
clc;

%% preparation 

% predefined parameters
freq = [8:1:15 8.2:1:15.2 8.4:1:15.4 8.6:1:15.6 8.8:1:15.8];
nCondition = 40;
channels = [1:9];% Pz, PO5, PO3, POz, PO4, PO6, O1, Oz, and O2
nChannel = length(channels);
rfs = 250;% sampling rate
dataLength = 6*rfs;% [-0.5 5.5s]
nBlock = 6;% six blocks 
latencyDelay = 0.5*rfs+35;% 140ms delay
stimTime = 5;% stimulation time (up to 5s)
N = round(stimTime*rfs);
t = [1:N]/rfs;

% predefined bandpass filter
fs = rfs/2;
Wp=[6/fs 90/fs];%6-90Hz
Ws=[4/fs 100/fs];%
[k,Wn]=cheb1ord(Wp,Ws,3,40);
[B,A] = cheby1(k,0.5,Wn);
freqz(B,A)

%% illustrate the PSD and SNR of SSVEP data

iTarget = 3;% 3:10Hz,8:15Hz
iChannel = 8;% 8:Oz
sTime = 5;% 5s stimulation

rawData=zeros(nChannel,dataLength,nCondition,nBlock);
% average across subjects
count=0;
for sub=[1:10]
   count=count+1;
   load(['S' int2str(sub) '.mat' ]);
   rawData = double(data(channels,1:dataLength,1:nCondition,1:nBlock));
   vepDataTmp = squeeze(mean(rawData(iChannel,latencyDelay+1:latencyDelay+sTime*rfs,iTarget,:),4));
   vepDataAll(count,:) = filtfilt(B,A,vepDataTmp);% filtering
   vepAmpAll(count,:) = abs(fft(vepDataAll(count,:)))*2/(length(vepDataAll(count,:)));%amplitude spectrum
end
vepData=mean(vepDataAll,1);
vepAmp=mean(vepAmpAll,1);

n=5;%width of background neighbours [-1 +1]Hz
vepSnr=zeros(length(vepAmp));
f1=6;%SNR starts from f1 Hz
for k=f1*sTime+1:length(vepAmp)-n
vepSnr(k)=20*log10(vepAmp(k)/mean(vepAmp([k-n:k-1 k+1:k+n])));%SNR
end


figure(2);
% temporal waveform
subplot(3,1,1);
plot([1:sTime*rfs]/rfs,vepData)
xlabel('Time(s)')
ylabel('Amplitude(uV)')

% amplitude spectrum
subplot(3,1,2);
plot(([1:sTime*rfs]-1)/sTime,vepAmp)
hold on
plot(freq(iTarget)*[1:7],vepAmp(freq(iTarget)*sTime*[1:7]+1),'ro')
xlim([f1 70])
xlabel('Frequency(Hz)')
ylabel('Amplitude(uV)')

% signal-to-noise ratio
subplot(3,1,3);
plot(([1:sTime*rfs]-1)/sTime,vepSnr)
hold on
plot(freq(iTarget)*[1:7],vepSnr(freq(iTarget)*sTime*[1:7]+1),'ro')
xlim([f1 70])
xlabel('Frequency(Hz)')
ylabel('SNR(dB)')
