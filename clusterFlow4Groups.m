function [idx,lengthPerCluster,curPoint,sizePerCluster]=clusterFlow4Groups(D,choice,clusterMax)
%map flow counters to groups
%: 1d, to k groups


%choice = 3; % 3, k-means, fixed cluster count;2, hierarchical cluster, max, 10 %1, affinity propagation, no fixed cluster count

% clusterMax = 10;
 
%%%%%%%%%
N = length(D);

if choice ==7
    
    %change
    
    %original, 90
   %compute the variance of each interval, and allocate the space by the interval.   
    % defaultOutlierThreshold = 90;
    defaultOutlierThreshold = 95;
    
%         cent=1;
%         sampledNum=ceil(cent*N);
%         tmp1=randperm(N);
%         %sample items
%         D = D(tmp1(1:sampledNum))';
%         clear tmp1;  
        
    %level
    pct90Val  = prctile(D,defaultOutlierThreshold);
    outlierVal = pct90Val;
  
    
    %level -1
    sampledItems = D(D<pct90Val);
   [IDX,C]=kmeans( sampledItems(:), clusterMax-1,...
                    'Distance','cityblock',...
                    'emptyaction','singleton',...
                    'Replicates',5,'display','off');
                
     clear sampledItems;           
    
     %sorted
     RowCentroids=sort(union(C,outlierVal));
    curPoint= sort(removeRedundantDat(RowCentroids));
    
    
    %find cluster index
    idx = zeros(N,1);
        for i=1:N
              Dtemp=abs(D(i)-curPoint);
              [BB,IX]=sort( Dtemp);       
               idx(i)=IX(1);            
        end;
   
  lengthPerCluster = zeros(length(curPoint),1);
  totalSum = 0;
    for clusterIdx=1:length(curPoint)   
        subsets = find(idx==clusterIdx);
         %fprintf(1,'%d %.2f %2f\n',length(subsets),std(D( subsets)),var(D( subsets)));
        if ~isempty(subsets)
            varsubSet = entropy4LSS(D( subsets));%need dimensionless, std, var(D( subsets));
            lengthPerCluster(clusterIdx) = varsubSet;
           
            totalSum = totalSum+varsubSet;
        end
    end
    %normalized
    for clusterIdx=1:length(curPoint)
         fprintf(1,'%d %d %2f %.2f\n',length(curPoint),clusterIdx,lengthPerCluster(clusterIdx),totalSum);
       lengthPerCluster(clusterIdx) =lengthPerCluster(clusterIdx) /totalSum;
    end
    
            if 1
        
    clusterIndexes = unique(idx)';    
    totalClusters=length(clusterIndexes);
    %cluster
    sizePerCluster = zeros(totalClusters,1);
    j=1;
    for i=clusterIndexes
          ii=find(idx==i);
          sizePerCluster(j) = length(ii)/length(D);
          j = j +1;
    end
            end
        
    
end

%header: <=2. tail, 90
if choice ==6
    
    firstRingThreshold = 2;
    
    defaultOutlierThreshold = 90;
    
%         cent=1;
%         sampledNum=ceil(cent*N);
%         tmp1=randperm(N);
%         %sample items
%         sampledItems0 = D(tmp1(1:sampledNum))';
%         clear tmp1;  
        
    %level
    pct90Val  = prctile(D,defaultOutlierThreshold);
    outlierVal = pct90Val;
    
    %level -1
    
    
    if(clusterMax>2)
        %level(2 -d-1)
        sampledItems0 = D(D<pct90Val);
        sampledItems=sampledItems0(sampledItems0>firstRingThreshold);
        clusterNum =clusterMax-2;
    else
        
        sampledItems=D;
         clusterNum =clusterMax;
    end
    
   [IDX,C]=kmeans( sampledItems(:), clusterNum ,...
                    'Distance','cityblock',...
                    'emptyaction','singleton',...
                    'Replicates',5,'display','off');
                
     clear sampledItems;     
     if(clusterMax>2)
     %add all cluster centers
        RowCentroids=sort(union(C,[outlierVal firstRingThreshold]));
        
     else
          RowCentroids=sort(C);
     end
    
    curPoint= removeRedundantDat(RowCentroids);
    idx = zeros(N,1);
        for i=1:N
              Dtemp=abs(D(i)-curPoint);
              [BB,IX]=sort( Dtemp);       
               idx(i)=IX(1);            
        end;

        if 1
        
    clusterIndexes = unique(idx)';    
    totalClusters=length(clusterIndexes);
    %cluster
    sizePerCluster = zeros(totalClusters,1);
    j=1;
    for i=clusterIndexes
          ii=find(idx==i);
          sizePerCluster(j) = length(ii)/length(D);
          j = j +1;
    end
        end
        
          lengthPerCluster = zeros(length(curPoint),1);
  totalSum = 0;
    for clusterIdx=1:length(curPoint)   
        subsets = find(idx==clusterIdx);
        if ~isempty(subsets)
            varsubSet = entropy4LSS(D( subsets));
            lengthPerCluster(clusterIdx) = varsubSet;
            totalSum = totalSum+varsubSet;
        end
    end
    %normalized
    for clusterIdx=1:length(curPoint)
       lengthPerCluster(clusterIdx) =lengthPerCluster(clusterIdx) /totalSum;
    end
    
    
end


%gaussian mixture model
if choice == 5
    
    %reset
    idx = [];
    %
    lengthPerCluster = zeros(clusterMax,1);
    
    %parameter
    defaultOutlierThreshold = 90;
    
        cent=1;
        sampledNum=ceil(cent*N);
        tmp1=randperm(N);
        %sample items
        D = D(tmp1(1:sampledNum))';
        clear tmp1;  
        
    %level
    pct90Val  = prctile(D,defaultOutlierThreshold);
    outlierVal = pct90Val;
    
    %length
    lengthPerCluster(clusterMax)= length(D(D>=pct90Val));
    
    %level -1
    sampledItems = D(D<pct90Val);
    %sort
    sampledItems = sort(sampledItems);
    sampledItems=sampledItems';
    
    %gaussian mixture model
    obj = fitgmdist(sampledItems,clusterMax-1,'RegularizationValue',0.1);
    %partition
    idxRandIndex = cluster(obj,sampledItems);
    
    fprintf(1,'gauss complete');
    
    %change idx, ascending order
    curPoint=zeros((clusterMax-1),1);
    for im=1:(clusterMax-1)
        result=sampledItems(idxRandIndex==im);
        %fprintf(1,'cluster: %d\n',length(result));
        if isempty(result)           
            continue;
        end
         curPoint(im) = mean(result);
         fprintf(1,'centroid: %d\n',curPoint(im));
    end
    
    fprintf(1,' curPoint complete');
    
    [curPoint,idxCluster] = sort(curPoint);
    
    fprintf(1,'resort centroid');
    %size of the corresponding point
    for im=1:(clusterMax-1)
              
        xx = idxCluster(im);       
        lengthPerCluster(xx) = length(sampledItems(idxRandIndex==im));
        fprintf(1,'index %d, total: %d, length: %d\n',im,length(idxCluster),lengthPerCluster(xx));
    end
    
    if 0    
    for ii=1:length(idxRandIndex)
        oldValue = idxRandIndex(ii);
        pos = find(idxCluster==oldValue);
        idxRandIndex(ii) = pos;
    end
    end
    %
    curPoint(clusterMax) = pct90Val;
     fprintf(1,'resort centroid');
    
    
    
end


if choice ==4 
       
    defaultOutlierThreshold = 90;
    
%         cent=1;
%         sampledNum=ceil(cent*N);
%         tmp1=randperm(N);
%         %sample items
%         D = D(tmp1(1:sampledNum))';
%         clear tmp1;  
        
    %level
    pct90Val  = prctile(D,defaultOutlierThreshold);
    outlierVal = pct90Val;
    %level -1
    sampledItems = D(D<pct90Val);
   [IDX,C]=kmeans( sampledItems(:), clusterMax-1,...
                    'Distance','cityblock',...
                    'emptyaction','singleton',...
                    'Replicates',5,'display','off');
                
     clear sampledItems;           
    RowCentroids=union(C,outlierVal);
    curPoint= removeRedundantDat(RowCentroids);
    idx = zeros(N,1);
        for i=1:N
              Dtemp=abs(D(i)-curPoint);
              [BB,IX]=sort( Dtemp);       
               idx(i)=IX(1);            
        end;

        if 1
        
    clusterIndexes = unique(idx)';    
    totalClusters=length(clusterIndexes);
    %cluster
    sizePerCluster = zeros(totalClusters,1);
    j=1;
    for i=clusterIndexes
          ii=find(idx==i);
          sizePerCluster(j) = length(ii)/length(D);
          j = j +1;
    end
        end
        
  lengthPerCluster = zeros(length(curPoint),1);
  totalSum = 0;
    for clusterIdx=1:length(curPoint)   
        subsets = find(idx==clusterIdx);
        if ~isempty(subsets)
            varsubSet = entropy4LSS(D( subsets));
            lengthPerCluster(clusterIdx) = varsubSet;
            totalSum = totalSum+varsubSet;
        end
    end
    %normalized
    for clusterIdx=1:length(curPoint)
       lengthPerCluster(clusterIdx) =lengthPerCluster(clusterIdx) /totalSum;
    end
    
    
end

if choice ==3

    
        cent=1;
        sampledNum=ceil(cent*N);
        tmp1=randperm(N);
        %sample items
        sampledItems = D(tmp1(1:sampledNum))';
        clear tmp1;  
        % k-means 
    [IDX,C]=kmeans( sampledItems(:), clusterMax,...
                    'Distance','cityblock',...
                    'emptyaction','singleton',...
                    'Replicates',5,'display','off');
                
     clear sampledItems;           
    RowCentroids=C';
    curPoint= removeRedundantDat(RowCentroids);
    
    idx = zeros(N,1);
        for i=1:N
              Dtemp=abs(D(i)-curPoint);
              [BB,IX]=sort( Dtemp);       
               idx(i)=IX(1);            
        end;

    totalClusters=length(unique(idx)');
    %cluster
    lengthPerCluster = zeros(totalClusters,1);
    j=1;
    for i=unique(idx)'
          ii=find(idx==i);
          lengthPerCluster(j) = length(ii);
          j = j +1;
    end
    
end

if choice ==2 
   
  
   M = pdist(D); 
   Z = linkage(M); 
   idx = cluster(Z,'maxclust',clusterMax); 
   totalClusters=length(unique(idx)');
    %cluster
    lengthPerCluster = zeros(totalClusters,1);
    j=1;
    for i=unique(idx)'
          ii=find(idx==i);
          lengthPerCluster(j) = length(ii);
          j = j +1;
    end


 end



if choice ==1
%affinity propagation: An algorithm that identifies exemplars among data points and forms 
%clusters of data points around these exemplars. It operates by simultaneously considering all data 
%point as potential exemplars and exchanging messages between data points until a good set of exemplars 
%and clusters emerges. (See BJ Frey and D Dueck, Science 315, Feb 16, 2007.
%

x=D; % Create N, 1-D data points
M=N*N-N; s=zeros(M,3); % Make ALL N^2-N similarities
j=1;
for i=1:N
  for k=[1:i-1,i+1:N]
    s(j,1)=i; s(j,2)=k; s(j,3)=-sum((x(i)-x(k)).^2);
    j=j+1;
  end;
end;
p=median(s(:,3)); % Set preference to median similarity
[idx,netsim,dpsim,expref]=apclusterSparse(s,p,'maxits',10,'plot');

fprintf('Number of clusters: %d\n',length(unique(idx)));
fprintf('Fitness (net similarity): %f\n',netsim);

totalClusters=length(unique(idx)');

lengthPerCluster = zeros(totalClusters,1);
j=1;
for i=unique(idx)'
      ii=find(idx==i);
      lengthPerCluster(j) = length(ii);
      j = j +1;
end

if 0
    
    figure; % Make a figures showing the data and the clusters
    for i=unique(idx)'
      ii=find(idx==i); h=plot(x(ii),x(ii),'o'); hold on;
      col=rand(1,3); set(h,'Color',col,'MarkerFaceColor',col);
      xi1=x(i,1)*ones(size(ii)); xi2=x(i)*ones(size(ii)); 
      line([x(ii),xi1]',[x(ii),xi2]','Color',col);
    end;

end

end
%repair NAN
lengthPerCluster = repairNAN(lengthPerCluster);
