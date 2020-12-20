function [precisionVal,recallVal]= F1Score4Vec(hhThreshold,estimatedFlowCounts, groundtruthFlowCounts)
%F1: 
%PR precision rate, ratio of true instances reported
%RR recall rate, ratio of reported true instances


%
N = length(groundtruthFlowCounts);

%index raw
[a,rawHHIndex] = sort(groundtruthFlowCounts,'descend');
%selected
x0= a>hhThreshold;

rawHHIndex0 = rawHHIndex(x0);

%index estimated
[a2,estHHIndex] = sort(estimatedFlowCounts,'descend');

x0= a2>hhThreshold;

estHHIndex0 = estHHIndex(x0);

precisionVal = length(intersect(rawHHIndex0,estHHIndex0))/length(estHHIndex0);

recallVal = length(intersect(rawHHIndex0,estHHIndex0))/length(rawHHIndex0);


