%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Hard output PSP
% Input:
%   yOut: received signal
%   stepSize: step size of LMS for channel update
%   constellation: modulate signal constellation
%   responseDuration: channel length
%   fIni: Initial channel coefficients
%   frequencyOffset1,frequencyOffset2: normalized frequency offsets of signals 1 and 2
%   CFOIndex: Initial phase index (to keep phase continuity for adjacent
%   batches)
% Output:
%   s1Decode, s2Decode: detected signals 1 and 2
%   fEstimate: estimated channel coefficients
%   eArray: distance between the received and estimated signals
% Edited by: Xiaobei
% 21/07/2021
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [s1Decode,s2Decode,fEstimate,metricMin]=psp2sources_fast(yOut,stepSize,constellation,responseDuration,fIni,frequencyOffset1,frequencyOffset2,CFOIndex)

M=length(constellation);

N=length(yOut);

numberState=M^(2*responseDuration);

% Initialization for viterbi algorithm
[stateMatrix,preStateMatrix]=iniViterbi2sources(M,responseDuration,constellation);

% Matrix to store path for each state
pathMatrix=zeros(numberState,N);

% Matrix to store metric for each state
metricMatrix=zeros(numberState,1);

eArray=zeros(1,N);

% Channel coefficients
fEstimate=fIni;
fEstimate1=fIni(1:2*responseDuration+1);
fEstimate2=fIni(2*responseDuration+2:end);

numberPrestate=M^2;

skMinus1=zeros(numberState*numberPrestate,responseDuration);
skMinus2=zeros(numberState*numberPrestate,responseDuration); 

CFO_phase1=exp(2*pi*1j*(CFOIndex+1:CFOIndex+N)*frequencyOffset1).'; 
CFO_phase2=exp(2*pi*1j*(CFOIndex+1:CFOIndex+N)*frequencyOffset2).'; 

seed1=reshape(repmat(constellation,M,1),1,M^2);
seed2=repmat(constellation,1,M);

for k=1:1:N
    skZero1=repmat(seed1,numberState,1);
    
    skZero2=repmat(seed2,numberState,1);
    
    skZero1=reshape(skZero1,numberState*numberPrestate,1);
    
    skZero2=reshape(skZero2,numberState*numberPrestate,1);
    
    skPlus1=repmat(stateMatrix(:,1:2:2*responseDuration-1),numberPrestate,1);
    
    skPlus2=repmat(stateMatrix(:,2:2:2*responseDuration),numberPrestate,1);
    
    stateMinus=preStateMatrix;
    for j=1:responseDuration 
        if k-j>0
            stateMinus=pathMatrix(stateMinus,k-j);
            skMinus1(:,responseDuration-j+1)=stateMatrix(stateMinus,1);
            skMinus2(:,responseDuration-j+1)=stateMatrix(stateMinus,2);
        else
            skMinus1(:,responseDuration-j+1)=0;
            skMinus2(:,responseDuration-j+1)=0;
        end
    end    
    
    % All possible combination of symbols from time index
    % k-responseDuration to k+responseDuration for signals 1 and 2
    sk1=[skMinus1,skZero1,skPlus1]*CFO_phase1(k);
    sk2=[skMinus2,skZero2,skPlus2]*CFO_phase2(k);
    
    % Estimated received signal for all possible combinations
    d=sk1*fEstimate1+sk2*fEstimate2;
    
    d=reshape(d,numberState,numberPrestate);
    
    d=d.';
          
    % Distance between received and estimated signals
    e=yOut(k)-d;
    
    % Branch metrics
    m=abs(e).^2;

    preMetric=metricMatrix(preStateMatrix).';    
    
    % Accumulated metrics
    gamma=preMetric+m;
    
    % Find the minimum accumulated metric and survivor path
    [metricMatrixUpdate,preStateIndex]=min(gamma,[],1); 
    
    i=1:numberState;
    
    index=(i-1)*size(e,1)+preStateIndex;    
    
    eSurvive=e(index).^2;
    
    preStateMatrixTmp=preStateMatrix.';
    
    preStateSurvive=preStateMatrixTmp(index);
    
    % Store survivor path
    pathMatrix(:,k)=preStateSurvive;
    
    % Update of the metric matrix
    metricMatrix=metricMatrixUpdate.';    
    
    % Update of the channel coefficients   
    tracebackNum=4;
    
    if k>2*responseDuration+tracebackNum-1
        [metricMin,stateNum]=min(metricMatrix);
        for j=k:-1:k-(2*responseDuration+tracebackNum-1)
            if stateNum>0
                s1Tmp(j-k+2*responseDuration+tracebackNum)=stateMatrix(stateNum,1);
                s2Tmp(j-k+2*responseDuration+tracebackNum)=stateMatrix(stateNum,2); 

                stateNum=pathMatrix(stateNum,j);
            else
                s1Tmp(j-k+2*responseDuration+tracebackNum)=0;
                s2Tmp(j-k+2*responseDuration+tracebackNum)=0; 
            end
        end
        
        s1Tmp(1:2*responseDuration+1)=s1Tmp(1:2*responseDuration+1)*CFO_phase1(k-tracebackNum);
        s2Tmp(1:2*responseDuration+1)=s2Tmp(1:2*responseDuration+1)*CFO_phase2(k-tracebackNum);
        
        x=[s1Tmp(1:2*responseDuration+1).';s2Tmp(1:2*responseDuration+1).'];
        dPrevious=fEstimate.'*x;         
        
        if responseDuration==1
            ePrevious=yOut(k-tracebackNum+1)-dPrevious;
        elseif responseDuration==2
            ePrevious=yOut(k-tracebackNum)-dPrevious;
        elseif responseDuration==3
            ePrevious=yOut(k-tracebackNum-1)-dPrevious;
        else
            error('responseDuration not valid');
        end         
           
        fEstimate=fEstimate+stepSize*ePrevious*conj(x); 
        fEstimate1=fEstimate(1:2*responseDuration+1);
        fEstimate2=fEstimate(2*responseDuration+2:end);   
    end      
    eArray(k)=min(eSurvive);    
end

s1Decode=zeros(N,1);
s2Decode=zeros(N,1);

% Trace back from the end to make decision
[metricMin,stateNum]=min(metricMatrix);
for k=N :-1:1
    if stateNum>0
        s1Decode(k)=stateMatrix(stateNum,1);
        s2Decode(k)=stateMatrix(stateNum,2);       
        
        stateNum=pathMatrix(stateNum,k);
    else
        s1Decode(k)=0;
        s2Decode(k)=0;
    end    
   
end



