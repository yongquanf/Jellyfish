function [FlowErr,F1Score,EntropyErr]=flowTelemetryAppsVaryThresholdJelly(D,sketchSize,choice,clusterMax, pSampled,useWeight,varyThreshold,DTrained, reverseIndex)
%telemetry

%1. flow distribution
[estimatedFlowCounts, groundtruthFlowCounts]=PredictFlowCount(D,sketchSize,choice,clusterMax, pSampled,useWeight,DTrained);

%reverse the elements
[estimatedFlowCounts]=reverseKeys(estimatedFlowCounts, reverseIndex);
[groundtruthFlowCounts]=reverseKeys(groundtruthFlowCounts, reverseIndex);

%relative error
%[RecCovs,REs,AEs]= distributionCovByRow(estimatedFlowCounts, groundtruthFlowCounts);
idx0 = groundtruthFlowCounts>0;
REs =  abs(groundtruthFlowCounts(idx0)-estimatedFlowCounts(idx0))./groundtruthFlowCounts(idx0);
clear idx0;


FlowErr=[mean(REs) median(REs) prctile(REs,90)];

clear REs AEs RecCovs;

%2. heavy hitter
%a flow whose byte count exceeds a threshold in an epoch
%estimate the flow count, and sort them, use F1 score to calculate the hh
%detect

%hh
if 0
percentHH = 0.001;
totalValues = sum(groundtruthFlowCounts);
hhThreshold = totalValues*percentHH ;
end
%10% is hh
hhThreshold = prctile(groundtruthFlowCounts,varyThreshold);

%compute the f1 score
%F1: 
%PR precision rate, ratio of true instances reported
%RR recall rate, ratio of reported true instances
[PR,RR]= F1Score4Vec(hhThreshold,estimatedFlowCounts, groundtruthFlowCounts);

%F1 score
F1Score = 2*PR*RR/(PR+RR);


%3. entropy
m = sum(groundtruthFlowCounts);

trueEntropy=0;
for i=1:length(groundtruthFlowCounts)
    tmp = groundtruthFlowCounts(i)/m;
    if tmp ==0 || isnan(tmp*log(tmp))
        continue;
    end
   trueEntropy = trueEntropy - tmp*log(tmp); 
end
m0 = sum(estimatedFlowCounts);
EstEntropy = 0;
for i=1:length(estimatedFlowCounts)
    tmp = estimatedFlowCounts(i)/m0;
    if tmp ==0 || isnan(tmp*log(tmp))
        continue;
    end
    EstEntropy = EstEntropy - tmp*log(tmp); 
end
 %entropy error
 EntropyErr = abs(EstEntropy - trueEntropy)/trueEntropy;





