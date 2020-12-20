function [errorDist,dists,idx]=distributionErrorPot(groundtruthFlowCounts,estimatedFlowCounts)
%divide to bins, for each bin, compute the mean, std

[~,idx] = sort(groundtruthFlowCounts);
groundtruthFlowCounts = groundtruthFlowCounts(idx);
estimatedFlowCounts = estimatedFlowCounts(idx);

dists=10:10:100;
curPoint = prctile(groundtruthFlowCounts,dists);
 m = length(curPoint);       
errorDist=zeros(m,2);   
for i=1:m
   point = curPoint(i);
   idx0 = groundtruthFlowCounts <= point;
   %filter
   if i>1
       leftPoint = curPoint(i-1);
       idx1=groundtruthFlowCounts> leftPoint;
       idx0 = idx0 & idx1;
       fprintf(1,'%.2f %.2f %d\n',point,curPoint(i-1),length(idx0(idx0>0)));
   else
       fprintf(1,'%.2f %d\n',point,length(idx0(idx0>0)));
   end
   
   fprintf(1,'#%d %d\n',length(groundtruthFlowCounts(idx0)),length(estimatedFlowCounts(idx0)));
   if isempty(groundtruthFlowCounts(idx0))
       continue;
   end
   %RE
   REs = abs(groundtruthFlowCounts(idx0)-estimatedFlowCounts(idx0))./groundtruthFlowCounts(idx0);
   
    errorDist(i,:)=[mean(REs) std(REs)];
    clear REs;
end


if 0
N = length(groundtruthFlowCounts);
%percentile map
    idx = zeros(N,1);
        for i=1:N
              %i
              Dtemp=abs(groundtruthFlowCounts(i)-curPoint);
              [BB,IX]=sort( Dtemp);   
              %
              sortPos=IX(1);   
              %pos
              %fprintf(1,'%d %d %d\n',i,sortPos,length(curPoint));
              
              if groundtruthFlowCounts > curPoint(sortPos)
                   sortPos=sortPos-1;
              end
              idx(i) = sortPos;
        end;
        
 m = length(curPoint);       
errorDist=zeros(m,2);        
tmp=cell(m,1);
for i=1:m
    tmp{i}=[];
end
    
for i=1:N
   groundTruthIndex = idx(i); 
   RE = abs(groundtruthFlowCounts(i)-estimatedFlowCounts(i))/groundtruthFlowCounts(i);
   %fprintf(1,'%d %.4f\n',groundTruthIndex,RE);
   array = tmp{groundTruthIndex};
   if isempty(array)
       array = RE;
       tmp{groundTruthIndex} = array;
   else
       array =[array RE];
       tmp{groundTruthIndex} = array;
   end
end

for i=1:m
    REs = tmp{i};
    errorDist(i,:)=[mean(REs) std(REs)];
end

clear tmp;
end

