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
weights = [1:7].^(-1.25)+0.25;
stimTime = 5;% stimulation time (up to 5s)
N = round(stimTime*rfs);
t = [1:N]/rfs;

% predefined bandpass filters for a filter bank
fs = rfs/2;
%h1
fls1(1) = [6]; 
fls2(1) = [4];
fhs1(1) = [90]; 
fhs2(1) = [100];
%h2
fls1(2) = [14]; 
fls2(2) = [10];
fhs1(2) = [90]; 
fhs2(2) = [100];
%h3
fls1(3) = [22]; 
fls2(3) = [16];
fhs1(3) = [90]; 
fhs2(3) = [100];
%h4
fls1(4) = [30]; 
fls2(4) = [24];
fhs1(4) = [90]; 
fhs2(4) = [100];
%h5
fls1(5) = [38]; 
fls2(5) = [32];
fhs1(5) = [90]; 
fhs2(5) = [100];
%h6
fls1(6) = [46]; 
fls2(6) = [40];
fhs1(6) = [90]; 
fhs2(6) = [100];
%h7
fls1(7) = [54]; 
fls2(7) = [48];
fhs1(7) = [90]; 
fhs2(7) = [100];
%Wn=0;
for nFB=1:7
Wp=[fls1(nFB)/fs fhs1(nFB)/fs];%
Ws=[fls2(nFB)/fs fhs2(nFB)/fs];%
[k,Wn]=cheb1ord(Wp,Ws,3,40);
[B{nFB},A{nFB}] = cheby1(k,0.5,Wn);
end

%% simulate online process
load(['S1.mat']);
rawData = double(data(channels,1:dataLength,1:nCondition,1:nBlock));

for ii = 1:nBlock% block 1-6

    % leave-one-out training data
    rawData1=rawData;
    rawData1(:,:,:,ii)=[];
    trainData(:,:,:)=mean(rawData1,4);
    for jj=1:nCondition
        for chan = 1:nChannel
        refData(chan,:,jj) = trainData(chan,1:N+latencyDelay,jj);
        refDatah1(chan,:,jj) = filtfilt(B{1},A{1},refData(chan,:,jj));
        refDatah2(chan,:,jj) = filtfilt(B{2},A{2},refData(chan,:,jj));
        refDatah3(chan,:,jj) = filtfilt(B{3},A{3},refData(chan,:,jj));
        refDatah4(chan,:,jj) = filtfilt(B{4},A{4},refData(chan,:,jj));
        refDatah5(chan,:,jj) = filtfilt(B{5},A{5},refData(chan,:,jj));
        refDatah6(chan,:,jj) = filtfilt(B{6},A{6},refData(chan,:,jj));
        refDatah7(chan,:,jj) = filtfilt(B{7},A{7},refData(chan,:,jj));
        end
    end
    
    for jj = 1:nCondition% target 1-40
        
           %tic
            % ----------------filter bank analysis-----------------------
            for chan = 1:nChannel
                testData = rawData(chan,1:N+latencyDelay,jj,ii);
                bpDatah1(chan,:) = filtfilt(B{1},A{1},testData);
                bpDatah2(chan,:) = filtfilt(B{2},A{2},testData);
                bpDatah3(chan,:) = filtfilt(B{3},A{3},testData);
                bpDatah4(chan,:) = filtfilt(B{4},A{4},testData);
                bpDatah5(chan,:) = filtfilt(B{5},A{5},testData);
                bpDatah6(chan,:) = filtfilt(B{6},A{6},testData);
                bpDatah7(chan,:) = filtfilt(B{7},A{7},testData);
            end
           
            % extract testdata epochs before CCA
            testDatah1 = (squeeze(bpDatah1(:,1+latencyDelay:N+latencyDelay)))';
            testDatah2 = (squeeze(bpDatah2(:,1+latencyDelay:N+latencyDelay)))';
            testDatah3 = (squeeze(bpDatah3(:,1+latencyDelay:N+latencyDelay)))';
            testDatah4 = (squeeze(bpDatah4(:,1+latencyDelay:N+latencyDelay)))';
            testDatah5 = (squeeze(bpDatah5(:,1+latencyDelay:N+latencyDelay)))';
            testDatah6 = (squeeze(bpDatah6(:,1+latencyDelay:N+latencyDelay)))';
            testDatah7 = (squeeze(bpDatah7(:,1+latencyDelay:N+latencyDelay)))';
            
            % IT-CCA for all sub-bands
            for kk = 1:nCondition
                refTmph1 = refDatah1(:,1+latencyDelay:N+latencyDelay,kk)';
                [SFA,SFB,D] = canoncorr(testDatah1,refTmph1);
                cc(1,kk)=corr2(testDatah1*SFA(:,1),refTmph1*SFA(:,1));
            
                refTmph2 = refDatah2(:,1+latencyDelay:N+latencyDelay,kk)';
                [SFA,SFB,D] = canoncorr(testDatah2,refTmph2);
                cc(2,kk)=corr2(testDatah2*SFA(:,1),refTmph2*SFA(:,1));

                refTmph3 = refDatah3(:,1+latencyDelay:N+latencyDelay,kk)';
                [SFA,SFB,D] = canoncorr(testDatah3,refTmph3);
                cc(3,kk)=corr2(testDatah3*SFA(:,1),refTmph3*SFA(:,1));

                refTmph4 = refDatah4(:,1+latencyDelay:N+latencyDelay,kk)';
                [SFA,SFB,D] = canoncorr(testDatah4,refTmph4);
                cc(4,kk)=corr2(testDatah4*SFA(:,1),refTmph4*SFA(:,1));

                refTmph5 = refDatah5(:,1+latencyDelay:N+latencyDelay,kk)';
                [SFA,SFB,D] = canoncorr(testDatah5,refTmph5);
                cc(5,kk)=corr2(testDatah5*SFA(:,1),refTmph5*SFA(:,1));

                refTmph6 = refDatah6(:,1+latencyDelay:N+latencyDelay,kk)';
                [SFA,SFB,D] = canoncorr(testDatah6,refTmph6);
                cc(6,kk)=corr2(testDatah6*SFA(:,1),refTmph6*SFA(:,1));

                refTmph7 = refDatah7(:,1+latencyDelay:N+latencyDelay,kk)';
                [SFA,SFB,D] = canoncorr(testDatah7,refTmph7);
                cc(7,kk)=corr2(testDatah7*SFA(:,1),refTmph7*SFA(:,1));
               
            end
             
            % weighted sum          
            rrr = weights*(sign(cc).*cc.^2);
            
            % target identification
            pLabels(ii,jj) = find(rrr==max(rrr));
        
            disp(  ['Block ' int2str(ii) ' - Target ' int2str(jj) ': estimated as ' int2str(pLabels(ii,jj)) ])
            %toc
        end
end

%% calculate accuracy

tLabels = repmat(1:nCondition,nBlock,1);
errors = pLabels-tLabels;
Accuracy = length(find(errors(:)==0))/(nCondition*nBlock)*100;
fprintf('Accuracy = %2.2f%% \n',Accuracy);
switchTime = 0.5;
Itr = itr(nCondition,Accuracy/100,stimTime+switchTime);
fprintf('ITR = %2.2f bpm \n\n',Itr);

