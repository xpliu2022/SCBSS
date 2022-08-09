%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initialize the Viterbi algorithm
% Input:
%   M: modulation order
%   responseDuration: channel length
%   constellation: modulate signal constellation
% Output:
%   stateMatrix: trellis states in Viterbi Algorithm
%   preStateMatrix: previous state for each state in stateMatrix
% Edited by: Xiaobei
% 21/07/2021
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [stateMatrix,preStateMatrix]=iniViterbi2sources(M,responseDuration,constellation)

numberState=M^(2*responseDuration);

stateMatrix=zeros(numberState,responseDuration);

%construct state matrix in Viterbi
for i=1:2*responseDuration
    seedMat=repmat(constellation,M^(i-1),1);
    seedVec=reshape(seedMat,M^i,1);
    stateColumn=repmat(seedVec,M^(2*responseDuration-i),1);
    stateMatrix(:,2*responseDuration-i+1)=stateColumn;
end

%previous states of each state
tempMat=repmat(1:numberState,M^2,1);
tempVec=reshape(tempMat,numberState*M^2,1);
preStateMatrix=reshape(tempVec,numberState,M^2);

end

    
        

        

