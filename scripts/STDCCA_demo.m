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

% reference signals for CCA
for ii = 1:length(freq);
    s1 = sin(2*pi*freq(ii)*t);
    s2 = cos(2*pi*freq(ii)*t);
    s3 = sin(2*pi*2*freq(ii)*t);
    s4 = cos(2*pi*2*freq(ii)*t);
    s5 = sin(2*pi*3*freq(ii)*t);
    s6 = cos(2*pi*3*freq(ii)*t);
    s7 = sin(2*pi*4*freq(ii)*t);
    s8 = cos(2*pi*4*freq(ii)*t);
    s9 = sin(2*pi*5*freq(ii)*t);
    s10 = cos(2*pi*5*freq(ii)*t);
    refData(:,:,ii) = cat(2,s1',s2',s3',s4',s5',s6',s7',s8',s9',s10');   
end
figure;
for i=1:10
subplot(10,1,i);plot(t,refData(:,i,3));
end
nHarmonics=5;

% predefined bandpass filter
fs = rfs/2;
Wp=[6/fs 90/fs];%6-90Hz
Ws=[4/fs 100/fs];%
[k,Wn]=cheb1ord(Wp,Ws,3,40);
[B,A] = cheby1(k,0.5,Wn);

%% simulate online process

load(['S1.mat']);
rawData = double(data(channels,1:dataLength,1:nCondition,1:nBlock));

for ii = 1:nBlock% block 1-6
 
    for jj = 1:nCondition% target 1-40
        
           %tic
            % ----------------filtering----------------------
            for chan = 1:nChannel
                testData = rawData(chan,1:N+latencyDelay,jj,ii);
                bpData(chan,:) = filtfilt(B,A,testData);
            end
           
            % extract testdata epochs before CCA
            testData = (squeeze(bpData(:,1+latencyDelay:N+latencyDelay)))';
            
            % CCA
            for kk = 1:nCondition
                [~,~,D] = canoncorr(testData, refData(:,1:2*nHarmonics,kk));
                cc(kk)=D(1);
             end
            
            % target identification
            pLabels(ii,jj) = find(cc==max(cc));
        
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

