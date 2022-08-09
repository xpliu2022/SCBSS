%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compare two symbol sequences and get the relative time and phase shift of the first to the second sequence 
% Input:
%   s1Decode: detected signal 1
%   s2Decode: detected signal 2
%   yOut: received signal
%   fIni: Initial channel coefficients
%   stepSize: step size of LMS for channel update
%   responseDuration: channel length
%   frequencyOffset1,frequencyOffset2: normalized frequency offsets of signals 1 and 2
%   CFOIndex: Initial phase index (to keep phase continuity for adjacent
%   batches
%   overlapLen: batch overlap length)
% Output:
%   fOut: estimated channel coefficients
% Edited by: Xiaobei
% 31/12/2021
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function fOut=fTrain(s1Decode,s2Decode,yOut,fIni,stepSize,responseDuration,frequencyOffset1,frequencyOffset2,CFOIndex,overlapLen)

% initialize the channel coefficients
fEstimate=fIni;
N=min(length(s1Decode),length(s2Decode));

%Determine the offset value to align the received signal yOut with estimated signal yReceive
if responseDuration==1
    offset=2;
else
    offset=3;
end

e=zeros(1,N-2*responseDuration);
CFO_phase1=exp(2*pi*1j*(CFOIndex+1:CFOIndex+N)*frequencyOffset1).'; 
CFO_phase2=exp(2*pi*1j*(CFOIndex+1:CFOIndex+N)*frequencyOffset2).'; 

% Use LMS for channel estimation
s1Decode(1:N)=s1Decode(1:N).*CFO_phase1;
s2Decode(1:N)=s2Decode(1:N).*CFO_phase2;
for i=1:N-overlapLen
    x=[s1Decode(i:i+2*responseDuration);s2Decode(i:i+2*responseDuration)];

    % estimate the received signal
    yReceive=fEstimate.'*x;

    % compute the error
    e(i)=yOut(i+offset)-yReceive;

    %update the estimated channel coefficients
    fEstimate=fEstimate+stepSize*e(i)*conj(x);
end 

f1=fEstimate(1:2*responseDuration+1);
f2=fEstimate(2*responseDuration+2:end);

% find the peak place of the channel coefficients
[f1_max,f1_max_index]=max(abs(f1));
[f2_max,f2_max_index]=max(abs(f2));

% shift the channel coefficients so that the peak of the channel coefficients will be in the center of the
% coefficient vector. For example, the estimated coefficients are [0.3 0.8
% 0.6 0.5 0.4], after shifting, it becomes [0.4 0.3 0.8 0.6 0.5]
f1=circshift(f1,responseDuration+1-f1_max_index,1);
f2=circshift(f2,responseDuration+1-f2_max_index,1);

% set the first and last channel coefficients to 0. This step is to make
% sure that the first and last channel coefficients won't be too large
f1(1)=0;
f1(end)=0;

f2(1)=0;
f2(end)=0;

fOut=[f1;f2];

end
    
    
    
