% This program is to simulate soft output per-survivor-processing (PSP)
% algorithm for co-channel signals detection. Inside hard output PSP,
% Viterbi algorithm (VA) is used.
% Edited by: Xiaobei
% 21/07/2021

clc;
clear all;
close all;

EbN0db=5:5:20;  % Eb/No

factor=1;

dataLen=1000*factor;    % data length

overlapLen=100*factor; % overlap length

tailLen=10*factor; % tail length

frontOffset=6;

blockLen=overlapLen+dataLen+tailLen; % block length

blockNum=25; % total block number

N=dataLen*(blockNum+1)+overlapLen+tailLen; % total block length

testNum=10; % test number

modType='QPSK'; % modulation type 

stepSize=0.01;  % stepSize in PSP

powerDiff=2; % power difference between the two signals, needs to be set to be > 0

tau = 3; % fractional delay (samples) between 2 signal

txOverSampling=10; % oversampling ratio

rxOverSampling=2; 

A1=1; % amplitude of signal 1
A2=sqrt(A1/(10^(powerDiff/10))); % amplitude of signal 2

% modulation order
if strcmp(modType,'BPSK')==1
    modOrder=2;
elseif strcmp(modType,'QPSK')==1
    modOrder=4;    
elseif strcmp(modType,'8PSK')==1
    modOrder=8;
else
    error('invalid modulation type');
end
bitPerSym=log2(modOrder); % number of bits per symbol

responseDuration=2; % channel length

numberState=(modOrder)^(2*responseDuration); %number of states in Viterbi algorithm

sigMeanPower=1;

% parameters for raised cosine filter
rrcSpan = 32; % span of symbol;
rolloff = 0.25; % Rolloff factor of filter
rrcfilter = rcosdesign(rolloff, rrcSpan, txOverSampling,'sqrt');  % Create a square root raised cosine filter.
filtOrder = rrcSpan;
filtDelay = filtOrder;

addFrequencyOffset=1; % 1: add frequency offset, 0: don't add

% frequency offsets for the 2 signals
if addFrequencyOffset==1
    frequencyOffset1=20;
    frequencyOffset2=50;
else
    frequencyOffset1=0;
    frequencyOffset2=0;
end

phi1 = 0;% initial phase for signal 1        
phi2 = 0;% initial phase for signal 2

symRate=10000; % symbol rate

%% start testing
errNum1_sum=zeros(1,length(EbN0db));
errNum2_sum=zeros(1,length(EbN0db));

for i=1:length(EbN0db)
     
    SNR=EbN0db(i)+10*log10(bitPerSym); %SNR of the signal 
          
    noiseDeviation = 1/(10^(SNR/20));     
    
    errNumSum1=0;
    errNumSum2=0;   
       
    tic
    
    parfor frameNum=1:testNum      
        % Generate source data          
        s1 = randi([0,1],bitPerSym*N,1); 
        s2 = randi([0,1],bitPerSym*N,1);  
       
        % Convert the bits in x into k-bit symbols for modulation.
        tx1 = bi2de(reshape(s1,bitPerSym,N).','left-msb');        
        tx2 = bi2de(reshape(s2,bitPerSym,N).','left-msb');     
        
        % Modulate the source signal        
        sModulate1 = pskmod(tx1,modOrder,0,'gray');
        sModulate2 = pskmod(tx2,modOrder,0,'gray');
        
        constellation = pskmod((0:modOrder-1),modOrder,0,'gray');          
               
         %% Add frequency offset
        t=0:1:N-1;
        t=t/symRate;
        
        CFO_phase1=exp(2*pi*1j*t*frequencyOffset1).';         
        CFO_phase2=exp(2*pi*1j*t*frequencyOffset2).';   

        % Add initial phase
        sModulateCFO1=sModulate1.*CFO_phase1*exp(1j*phi1);
        sModulateCFO2=sModulate2.*CFO_phase2*exp(1j*phi2);
        
        % Pass signal through rrc filter
        sFlt1 = A1*upfirdn(sModulateCFO1, rrcfilter, txOverSampling);        
        sFlt2 = A2*upfirdn(sModulateCFO2, rrcfilter, txOverSampling);
        
        % Introduce non-integer OSR
        sFlt1 = resample(sFlt1,9999,10000);
        sFlt2 = resample(sFlt2,10001,10000);
        
        % Pass signal through channel
        s_Ch1=[sFlt1;zeros(tau*txOverSampling,1)];
        s_Ch2=[zeros(tau*txOverSampling,1);sFlt2];     
        
%         s_Ch1=[sFlt1;zeros(tau,1)];
%         s_Ch2=[zeros(tau,1);sFlt2];     
               
        s_length=min(length(s_Ch1),length(s_Ch2));
        
        awgnDefault = sqrt(sigMeanPower)/sqrt(2)*(randn(s_length,1)  + sqrt(-1)*randn(s_length,1));
    
        noise = noiseDeviation*awgnDefault;
        
        % Received signal
        yrx = s_Ch1(1:s_length) + s_Ch2(1:s_length) + noise;       
        
        % Pass received signal through rrc filter
        yrx_filt = filter(rrcfilter,1,yrx);
        
        yrx_filt = yrx_filt((filtDelay-frontOffset)*txOverSampling+1:end);

        % Downsample the received signal
        yrx1 = yrx_filt(3:txOverSampling:end);
        
        % Remove the "0s" before signal arrival
        yrxDS=yrx1(filtDelay+1:end);  
                           
%% channel estimation and equalization        
 
        % Initialize the channel coefficients        
        fIni=zeros(2*(2*responseDuration+1),1);
        fIni(responseDuration+1)=1;%A1*exp(1j*phi1);
        fIni(3*responseDuration+2)=1;%A2*exp(1j*phi2);  
        
        fIni_data=repmat(fIni,1,rxOverSampling);
        
        eqSig1=zeros(blockNum*dataLen+overlapLen,1);
        eqSig2=zeros(blockNum*dataLen+overlapLen,1);        
            
        last_index1=0;
        last_index2=0;             
            
        % Start equalization batch-by-batch
        for j=1:blockNum
            
            startPoint=(j-1)*dataLen;
            
            yrx_resample=resample(yrx_filt,rxOverSampling,10);   
            
            yOut=zeros(blockLen,rxOverSampling);
            s1Decode=zeros(blockLen,rxOverSampling);
            s2Decode=zeros(blockLen,rxOverSampling);
            fEstimate=zeros(2*(2*responseDuration+1),rxOverSampling);
            m=zeros(1,rxOverSampling);
           
            for index=1:rxOverSampling
                yOut(:,index)=yrx_resample(startPoint*rxOverSampling+index:rxOverSampling:(startPoint+blockLen)*rxOverSampling);            
                                         
                % Equalization by PSP
                [s1Decode(:,index),s2Decode(:,index),fEstimate(:,index),m(index)]=psp2sources_fast(yOut(:,index),stepSize,constellation,responseDuration,fIni_data(:,index),frequencyOffset1/symRate,frequencyOffset2/symRate,startPoint); 

            end                  

            [min_m,min_m_index]=min(m);                

            s1Decode=s1Decode(:,min_m_index);  
            s2Decode=s2Decode(:,min_m_index);                       

            if j==1                               
                eqSig1(startPoint+1:startPoint+dataLen+overlapLen)=s1Decode(1:dataLen+overlapLen);
                eqSig2(startPoint+1:startPoint+dataLen+overlapLen)=s2Decode(1:dataLen+overlapLen);              

                last_index1=startPoint+dataLen+overlapLen;                
                last_index2=startPoint+dataLen+overlapLen;                                    
            else                
                [shiftNum1,phase1,peak1]=compareSequence(s1Decode(1:overlapLen),eqSig1(last_index1-overlapLen+1:last_index1),20,overlapLen-20,18);
                [shiftNum2,phase2,peak2]=compareSequence(s2Decode(1:overlapLen),eqSig2(last_index2-overlapLen+1:last_index2),20,overlapLen-20,18);           

                eqSig1(last_index1+1:last_index1+dataLen+1)=s1Decode(overlapLen+1+shiftNum1:overlapLen+dataLen+1+shiftNum1)*exp(1j*phase1);
                eqSig2(last_index2+1:last_index2+dataLen+1)=s2Decode(overlapLen+1+shiftNum2:overlapLen+dataLen+1+shiftNum2)*exp(1j*phase2); 


                 % to insure the shiftNum is not too small nor too big
                if shiftNum1<-5
                    last_index1=last_index1+dataLen+1;   
                elseif shiftNum1>5
                    last_index1=last_index1+dataLen-1;
                else
                    last_index1=last_index1+dataLen;
                end

                if shiftNum2<-5
                    last_index2=last_index2+dataLen+1; 
                elseif shiftNum2>5
                    last_index2=last_index2+dataLen-1;
                else
                    last_index2=last_index2+dataLen;
                end        
            end

            i
            j

            for index=1:rxOverSampling
               fIni_data(:,index)=fTrain(s1Decode,s2Decode,yOut(:,index),fIni,stepSize,responseDuration,frequencyOffset1/symRate,frequencyOffset2/symRate,startPoint,overlapLen);
            end
        end       
        
        % Time and phase synchronization (for BER calculation purpose,
        % should not be included in real application)
        [shiftNum1,phases1,peaks11]=compareSequence(eqSig1(1:dataLen*3),sModulate1(1:dataLen*3),overlapLen,dataLen*3-30,30);
        [shiftNum2,phases2,peaks22]=compareSequence(eqSig2(1:dataLen*3),sModulate2(1:dataLen*3),overlapLen,dataLen*3-30,30);
        
         % Keep the detected signal after synchronization
        eqSigData1=eqSig1(overlapLen+shiftNum1+1:overlapLen+shiftNum1+(blockNum-1)*dataLen); 
        eqSigData2=eqSig2(overlapLen+shiftNum2+1:overlapLen+shiftNum2+(blockNum-1)*dataLen);                
      
        % Demodulation
        demodData1=pskdemod(eqSigData1,modOrder,0,'gray').';
        demodData2=pskdemod(eqSigData2,modOrder,0,'gray').';
       
        demodBit1=reshape(de2bi(demodData1,'left-msb')',bitPerSym*(dataLen*(blockNum-1)),1);  
        demodBit2=reshape(de2bi(demodData2,'left-msb')',bitPerSym*(dataLen*(blockNum-1)),1);       
        
        % Calculate the number of errors before decoding
        [errNum1,errRate1,errPlace1]=biterr(demodBit1,s1(bitPerSym*overlapLen+1:bitPerSym*(overlapLen+dataLen*(blockNum-1))));
        [errNum2,errRate2,errPlace2]=biterr(demodBit2,s2(bitPerSym*overlapLen+1:bitPerSym*(overlapLen+dataLen*(blockNum-1))));        
        
%         figure(2)
%         plot(errPlace1);
        
%         save('errorPlace2.mat','errPlace1');
       
        errNumSum1=errNumSum1+errNum1;
        errNumSum2=errNumSum2+errNum2;
       
              
    end
    toc
    
    errNum1_sum(i) = errNumSum1;
    errNum2_sum(i) = errNumSum2; 
end

ber1=zeros(1,length(EbN0db));
ber2=zeros(1,length(EbN0db));

for i=1:length(EbN0db)    
    ber1(i)=errNum1_sum(i)/(testNum*bitPerSym*dataLen*(blockNum-1));
    ber2(i)=errNum2_sum(i)/(testNum*bitPerSym*dataLen*(blockNum-1));   
end

figure(1)
semilogy(EbN0db,ber1,'-^');
hold on
semilogy(EbN0db,ber2,'-o');
grid on
xlabel('Eb/No(dB)');
ylabel('BER');
legend('strong signal','weak signal','location','northeast');
hold off

