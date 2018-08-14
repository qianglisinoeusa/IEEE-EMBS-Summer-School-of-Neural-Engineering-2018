clear all;
clc;

%% preparation 

% predefined parameters
freq = [8:1:15 8.2:1:15.2 8.4:1:15.4 8.6:1:15.6 8.8:1:15.8];
nCondition = 40;
channels = [1:9];% Pz, PO5, PO3, POz, PO4, PO6, O1, Oz, and O2
nChannel = length(channels);
iChannel = [8];% Oz for FFT
rfs = 250;% sampling rate
dataLength = 6*rfs;% [-0.5 5.5s]
nBlock = 6;% six blocks 
latencyDelay = 0.5*rfs+35;% 140ms delay
stimTime = 5;% stimulation time (up to 5s)
N = round(stimTime*rfs);
t = [1:N]/rfs;

%% simulate online process

load(['S1.mat']);
rawData = double(data(channels,1:dataLength,1:nCondition,1:nBlock));

for ii = 1:nBlock% block 1-6
 
    for jj = 1:nCondition% target 1-40
        
           %tic
           
           % extract testdata epochs before FFT
            testData = (squeeze(rawData(iChannel,1+latencyDelay:N+latencyDelay,jj,ii)))';
            
            % FFT
            vepFFT=abs(fft(testData,5*rfs));%5*rfs-point FFT, zero-padding
            for kk = 1:nCondition
                amp(kk)=mean(vepFFT(freq(kk)*5*[1:5]+1));% mean of fundamental and harmonics
            end
           
            % target identification
            pLabels(ii,jj) = find(amp==max(amp));
        
            disp(  ['Block ' int2str(ii) ' - Target ' int2str(jj) ': estimated as ' int2str(pLabels(ii,jj)) ])
            %toc
        end
end

%% calculate accuracy

tLabels = repmat(1:nCondition,6,1);
errors = pLabels-tLabels;
Accuracy = length(find(errors(:)==0))/(nCondition*nBlock)*100;
fprintf('Accuracy = %2.2f%% \n',Accuracy);
switchTime = 0.5;
Itr = itr(nCondition,Accuracy/100,stimTime+switchTime);
fprintf('ITR = %2.2f bpm \n\n',Itr);

