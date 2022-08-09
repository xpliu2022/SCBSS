%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compare two symbol sequences and get the relative time and phase shift of the first to the second sequence 
% Input:
%   s1: first symbol sequence
%   s2: second symbol sequence
%   startPoint: start point of the sequence to compare
%   endPoint: end point of the sequence to compare
%   maxShift: maximum number of shift symbols to test
% Output:
%   shiftNum: relative time (symbols) shift of s1 compared to s2
%   phase: relative phase shift of s1 compared to s2
%   peak: maximum correlation between s1 and s2
% Edited by: Xiaobei
% 21/07/2021
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [shiftNum,phase,peak]=compareSequence(s1,s2,startPoint,endPoint,maxShift)

if length(s1)<endPoint
    error('s1 must be longer than endPoint');
end

N=length(s2);
if startPoint-floor(maxShift/2)<0||endPoint+floor(maxShift/2)>N  
    error('wrong setting');    
end

corr=zeros(maxShift,1);
phase_corr=zeros(maxShift,1);

% Calculate the correlation between s1 and s2
for i=1:1:maxShift
    temp=s1(startPoint+i+1-floor(maxShift/2):endPoint+i-floor(maxShift/2)).*conj(s2(startPoint+1:endPoint));
    corr(i)=abs(sum(temp));
    phase_corr(i)=angle(mean(temp));
end

% Find the peak of the correlation
[peak,peak_index]=max(corr);

% Find the relative time and phase shift
shiftNum=peak_index-floor(maxShift/2);
phase=-phase_corr(peak_index);



