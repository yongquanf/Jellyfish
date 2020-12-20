function [estimatedFlowCounts, groundtruthFlowCounts]=PredictFlowCount(D,sketchSize,choice,clusterMax, pSampled,useWeight,DTrained)
%choice: 0, rps, 1, cs, 2, cm

len0 = length(D);

%StringFormat=string('PredictFlowCoun_%d_Cluster_%d_UseWt_%d');

fnameA = sprintf('PredictFlowCoun_%d_Cluster_%d_UseWt_%d',choice,clusterMax,useWeight);

ExpLog = fopen(fnameA,'a');

%StringFormat=string();
fname0 = sprintf('PredictFlowCoun_%d_sketchSiz_%d_Dlen_%d_Cluster_%d_UseWt_%d',choice,sketchSize,len0,clusterMax,useWeight);
fprintf(ExpLog,'%s\n',fname0);

    N = length(D);
    d = sketchSize;
    %raw
    groundtruthFlowCounts = D;

 
        if choice ==-8
        
        %config
       %choiceCluster = 6; %, header, + tail ring, 
       %choiceCluster = 4; %k-means
       choiceCluster = 7; % assign, ratio of variance
       
         %cluster
        [idx,lengthPerCluster,curPoint,sizePerCluster]=clusterFlow4Groups(DTrained,choiceCluster,clusterMax);
        
        
        
        estimatedFlowCounts=zeros(N,1);

    %relative points
        sumCurPoints = sum(curPoint);
        curPointRelative = curPoint./sumCurPoints;
         %total eval
            totalVal = 0;
            for iL=1:length(lengthPerCluster)
               if useWeight ==0
                   totalVal= totalVal+lengthPerCluster(iL)*curPointRelative(iL)*sizePerCluster(iL);
                else if useWeight==1 
                    totalVal= totalVal+lengthPerCluster(iL)*curPointRelative(iL);
                    else if useWeight ==2
                            totalVal= totalVal+lengthPerCluster(iL)*sizePerCluster(iL);
                    
                            else if useWeight ==3
                                    totalVal= totalVal+lengthPerCluster(iL);
                                    else if useWeight ==4
                            
                                            totalVal= totalVal+ sizePerCluster(iL);
                                        else if useWeight ==5
                                            totalVal= totalVal+ curPointRelative(iL);
                                            else 
                                                totalVal= totalVal+ 1;
                                        end
                                end
                            end
                    end
               end
            end
        end
        
        %construct by flow distributions,
        %unique, also produce sorted value
        for i=unique(idx)'
            %index of cluster i
            ii=find(idx==i);
            %build sketch
            if 0
                sketchLetSize = max(ceil((length(ii)/sum(lengthPerCluster))*sketchSize),2);
            end
            

            %assign by variance
           
            if useWeight ==0
                normalizedVariance = lengthPerCluster(i)*curPointRelative(i)*sizePerCluster(i)/totalVal;
            else if useWeight==1
                normalizedVariance = (lengthPerCluster(i)*curPointRelative(i))/totalVal;
            else if useWeight==2
                    normalizedVariance = (lengthPerCluster(i)*sizePerCluster(i))/totalVal;
                else if useWeight==3
                    normalizedVariance = lengthPerCluster(i)/totalVal;
                    else if useWeight==4
                        normalizedVariance = sizePerCluster(i)/totalVal;
                     else if useWeight ==5
                        normalizedVariance = sizePerCluster(i)/totalVal;
                     else
                       normalizedVariance =  1/totalVal;
                    end
                end
                end
                    
                end
                end
            end
            
            %normalizedVariance? ratio
            
            if normalizedVariance<=0
                sketchLetSize =1;
            else
                sketchLetSize = ceil(normalizedVariance*sketchSize);
            end
            fprintf(1,'sketchlet: %.6f %.6f %.6f  %.6f %.6f %.6f \n',lengthPerCluster(i),curPointRelative(i),sizePerCluster(i),normalizedVariance,sketchSize,sketchLetSize);
            
            fprintf(ExpLog,'sketchlet: %.6f %.6f %.6f  %.6f %.6f %.6f \n',lengthPerCluster(i),curPointRelative(i),sizePerCluster(i),normalizedVariance,sketchSize,sketchLetSize);

            %estimatedSize
            [estimatedFlowCounts0]=PredictFlowCount(D(ii),sketchLetSize,1,clusterMax, pSampled,useWeight,DTrained);
            %fill
           estimatedFlowCounts(ii) = estimatedFlowCounts0;
        end
        
        
        end 
    
        
   
        if choice ==-7
        
        %config
       %choiceCluster = 6; %, header, + tail ring, 
       %choiceCluster = 4; %k-means
       choiceCluster = 7; % assign, ratio of variance
       
         %cluster
        [idx,lengthPerCluster,curPoint,sizePerCluster]=clusterFlow4Groups(D,choiceCluster,clusterMax);
        estimatedFlowCounts=zeros(N,1);
        
        % useWeight = 2;
         
         %total eval
         %relative points
    %relative points
        sumCurPoints = sum(curPoint);
        curPointRelative = curPoint./sumCurPoints;
         %total eval
            totalVal = 0;
            for iL=1:length(lengthPerCluster)
               if useWeight ==0
                   totalVal= totalVal+lengthPerCluster(iL)*curPointRelative(iL)*sizePerCluster(iL);
                else if useWeight==1 
                    totalVal= totalVal+lengthPerCluster(iL)*curPointRelative(iL);
                    else if useWeight ==2
                            totalVal= totalVal+lengthPerCluster(iL)*sizePerCluster(iL);
                    
                            else if useWeight ==3
                                    totalVal= totalVal+lengthPerCluster(iL);
                                    else if useWeight ==4
                            
                                            totalVal= totalVal+ sizePerCluster(iL);
                                        else if useWeight ==5
                                            totalVal= totalVal+ curPointRelative(iL);
                                            else 
                                                totalVal= totalVal+ 1;
                                        end
                                end
                            end
                    end
               end
            end
        end
        
        %construct by flow distributions,
        %unique, also produce sorted value
        for i=unique(idx)'
            %index of cluster i
            ii=find(idx==i);
            %build sketch
            if 0
                sketchLetSize = max(ceil((length(ii)/sum(lengthPerCluster))*sketchSize),2);
            end
            

            %assign by variance
           
                   if useWeight ==0
                normalizedVariance = lengthPerCluster(i)*curPointRelative(i)*sizePerCluster(i)/totalVal;
            else if useWeight==1
                normalizedVariance = (lengthPerCluster(i)*curPointRelative(i))/totalVal;
            else if useWeight==2
                    normalizedVariance = (lengthPerCluster(i)*sizePerCluster(i))/totalVal;
                else if useWeight==3
                    normalizedVariance = lengthPerCluster(i)/totalVal;
                    else if useWeight==4
                        normalizedVariance = sizePerCluster(i)/totalVal;
                     else if useWeight ==5
                        normalizedVariance = sizePerCluster(i)/totalVal;
                     else
                       normalizedVariance =  1/totalVal;
                    end
                end
                end
                    
                end
                end
            end
            
            %normalizedVariance? ratio
            
            if normalizedVariance<=0
                sketchLetSize =1;
            else
                sketchLetSize = ceil(normalizedVariance*sketchSize);
            end
            %fprintf(1,'sketchlet: %.2f %.2f  %.2f %.2f %.2f \n',lengthPerCluster(i),sizePerCluster(i),normalizedVariance,sketchSize,sketchLetSize);
            fprintf(1,'sketchlet: %.6f %.6f %.6f  %.6f %.6f %.6f \n',lengthPerCluster(i),curPointRelative(i),sizePerCluster(i),normalizedVariance,sketchSize,sketchLetSize);
           fprintf(ExpLog,'sketchlet: %.6f %.6f %.6f  %.6f %.6f %.6f \n',lengthPerCluster(i),curPointRelative(i),sizePerCluster(i),normalizedVariance,sketchSize,sketchLetSize);
          
            %estimatedSize
            [estimatedFlowCounts0]=PredictFlowCount(D(ii),sketchLetSize,2,clusterMax, pSampled,useWeight,DTrained);
            %fill
           estimatedFlowCounts(ii) = estimatedFlowCounts0;
        end
        
        
        end 
    
    
    
           %outlier+k-means-cluster + LCS + vary sampled percent
    if choice ==-6
        
        %config
        %clusterMax = 4;
        %choiceCluster = 4; %k-means
         %cluster
        % [idx,lengthPerCluster,curPoint,sizePerCluster]= clusterFlow4GroupsSampleRate(D, clusterMax, pSampled);
        %trained set
        [idx,lengthPerCluster,curPoint,sizePerCluster]= clusterFlow4GroupsVSampleSet(D,DTrained,clusterMax, pSampled);
        
         estimatedFlowCounts=zeros(N,1);
        %construct by flow distributions,
         %useWeight = 2;
         
         %total eval
     %relative points
        sumCurPoints = sum(curPoint);
        curPointRelative = curPoint./sumCurPoints;
         %total eval
            totalVal = 0;
            for iL=1:length(lengthPerCluster)
               if useWeight ==0
                   totalVal= totalVal+lengthPerCluster(iL)*curPointRelative(iL)*sizePerCluster(iL);
                else if useWeight==1 
                    totalVal= totalVal+lengthPerCluster(iL)*curPointRelative(iL);
                    else if useWeight ==2
                            totalVal= totalVal+lengthPerCluster(iL)*sizePerCluster(iL);
                    
                            else if useWeight ==3
                                    totalVal= totalVal+lengthPerCluster(iL);
                                    else if useWeight ==4
                            
                                            totalVal= totalVal+ sizePerCluster(iL);
                                        else if useWeight ==5
                                            totalVal= totalVal+ curPointRelative(iL);
                                            else 
                                                totalVal= totalVal+ 1;
                                        end
                                end
                            end
                    end
               end
            end
        end
        
        %construct by flow distributions,
        %unique, also produce sorted value
        for i=unique(idx)'
            %index of cluster i
            ii=find(idx==i);
            %build sketch
            if 0
                sketchLetSize = max(ceil((length(ii)/sum(lengthPerCluster))*sketchSize),2);
            end
            

            %assign by variance
                       if useWeight ==0
                normalizedVariance = lengthPerCluster(i)*curPointRelative(i)*sizePerCluster(i)/totalVal;
            else if useWeight==1
                normalizedVariance = (lengthPerCluster(i)*curPointRelative(i))/totalVal;
            else if useWeight==2
                    normalizedVariance = (lengthPerCluster(i)*sizePerCluster(i))/totalVal;
                else if useWeight==3
                    normalizedVariance = lengthPerCluster(i)/totalVal;
                    else if useWeight==4
                        normalizedVariance = sizePerCluster(i)/totalVal;
                     else if useWeight ==5
                        normalizedVariance = sizePerCluster(i)/totalVal;
                     else
                       normalizedVariance =  1/totalVal;
                    end
                end
                end
                    
                end
                end
            end
            
            %normalizedVariance? ratio
            
            if normalizedVariance<=0
                sketchLetSize =1;
            else
                sketchLetSize = ceil(normalizedVariance*sketchSize);
            end
            %fprintf(1,'sketchlet: %.2f %.2f  %.2f %.2f %.2f \n',lengthPerCluster(i),sizePerCluster(i),normalizedVariance,sketchSize,sketchLetSize);
            fprintf(ExpLog,'sketchlet: %.6f %.6f %.6f  %.6f %.6f %.6f \n',lengthPerCluster(i),curPointRelative(i),sizePerCluster(i),normalizedVariance,sketchSize,sketchLetSize);
          
            %estimatedSize
            [estimatedFlowCounts0]=PredictFlowCount(D(ii),sketchLetSize,0,clusterMax, pSampled,useWeight,DTrained);
            %fill
           estimatedFlowCounts(ii) = estimatedFlowCounts0;
        end
        
        if 0
        estimatedFlowCounts=zeros(N,1);
        %construct by flow distributions
        for i=unique(idx)'
            %index of cluster i
            ii=find(idx==i);
            %build sketch
            sketchLetSize = max(ceil((length(ii)/sum(lengthPerCluster))*sketchSize),2);
            %estimatedSize
            [estimatedFlowCounts0]=PredictFlowCount(D(ii),sketchLetSize,0);
            %fill
           estimatedFlowCounts(ii) = estimatedFlowCounts0;
        end
        end
        
        
    end 
    
    
       %outlier+k-means-cluster + rpsor
    if choice ==-5
        
        %config
       %choiceCluster = 6; %, header, + tail ring, 
       %choiceCluster = 4; %k-means
       choiceCluster = 7; % assign, ratio of variance
       
         %cluster, old
        %[idx,lengthPerCluster,curPoint,sizePerCluster]=clusterFlow4Groups(D,choiceCluster,clusterMax);
        %new, separate training set
        [idx,lengthPerCluster,curPoint,sizePerCluster]=clusterFlow4GroupsVSampleSet(D,DTrained,choiceCluster,clusterMax,ExpLog);
        
        estimatedFlowCounts=zeros(N,1);
        
        % useWeight = 2;
         
        %relative points
      %relative points
      %relative pos
      if 0
        sumCurPoints = sum(curPoint);
        curPointRelative = curPoint./sumCurPoints;
        
      end
      %direct
      curPointRelative = curPoint;
      
         %total eval
            totalVal = 0;
            for iL=1:length(lengthPerCluster)
               if useWeight ==0
                   totalVal= totalVal+lengthPerCluster(iL)*curPointRelative(iL)*sizePerCluster(iL);
                else if useWeight==1 
                        val = lengthPerCluster(iL)*curPointRelative(iL); 
                        if ~isnan(val) && val~=0
                            totalVal= totalVal+lengthPerCluster(iL)*curPointRelative(iL); 
                        end
                    else if useWeight ==2
                            totalVal= totalVal+lengthPerCluster(iL)*sizePerCluster(iL);
                    
                            else if useWeight ==3
                                    totalVal= totalVal+lengthPerCluster(iL);
                                    else if useWeight ==4
                            
                                            totalVal= totalVal+ sizePerCluster(iL);
                                        else if useWeight ==5
                                            totalVal= totalVal+ curPointRelative(iL);
                                            else 
                                                totalVal= totalVal+ 1;
                                        end
                                end
                            end
                    end
               end
            end
            end %end
        
        %construct by flow distributions,
        %unique, also produce sorted value
        for i=unique(idx)'
            %index of cluster i
            ii=find(idx==i);
            
            if isempty(ii)
                fprintf(1,'empty set skip: %d\n',i);
                continue;
            end
            %build sketch
            if 0
                sketchLetSize = max(ceil((length(ii)/sum(lengthPerCluster))*sketchSize),2);
            end
            

            %assign by variance
           
            if useWeight ==0
                normalizedVariance = lengthPerCluster(i)*curPointRelative(i)*sizePerCluster(i)/totalVal;
            else if useWeight==1
                normalizedVariance = (lengthPerCluster(i)*curPointRelative(i))/totalVal;
            else if useWeight==2
                    normalizedVariance = (lengthPerCluster(i)*sizePerCluster(i))/totalVal;
                else if useWeight==3
                    normalizedVariance = lengthPerCluster(i)/totalVal;
                    else if useWeight==4
                        normalizedVariance = sizePerCluster(i)/totalVal;
                     else if useWeight ==5
                        normalizedVariance = sizePerCluster(i)/totalVal;
                     else
                       normalizedVariance =  1/totalVal;
                    end
                end
                end
                    
                end
                end
            end
            
            %normalizedVariance? ratio
            
            if normalizedVariance<=0 || isnan(normalizedVariance)
                sketchLetSize =1; %old 1
            else
                sketchLetSize = max(1,ceil(normalizedVariance*sketchSize));
            end
            %fprintf(1,'sketchlet: %.2f %.2f %.2f  %.2f %.2f %.2f \n',lengthPerCluster(i),sizePerCluster(i),curPointRelative(i),normalizedVariance,sketchSize,sketchLetSize);
             fprintf(1,'sketchlet: %.6f %.6f %.6f  %.6f %.6f %.6f \n',lengthPerCluster(i),curPointRelative(i),sizePerCluster(i),normalizedVariance,sketchSize,sketchLetSize);
            fprintf(ExpLog,'sketchlet: %.6f %.6f %.6f  %.6f %.6f %.6f \n',lengthPerCluster(i),curPointRelative(i),sizePerCluster(i),normalizedVariance,sketchSize,sketchLetSize);
          
            %estimatedSize
            [estimatedFlowCounts0]=PredictFlowCount(D(ii),sketchLetSize,0,clusterMax, pSampled,useWeight,DTrained);
            %fill
           estimatedFlowCounts(ii) = estimatedFlowCounts0;
        end
        
        
    end 
    
     %cluster + cm
        if choice ==-4
        
        %config
        %clusterMax = 4;
        choiceCluster = 3; %k-means
         %cluster
        [idx,lengthPerCluster]=clusterFlow4Groups(D,choiceCluster,clusterMax);
        estimatedFlowCounts=zeros(N,1);
        %construct by flow distributions
        for i=unique(idx)'
            %index of cluster i
            ii=find(idx==i);
            %build sketch
            sketchLetSize = max(ceil((length(ii)/sum(lengthPerCluster))*sketchSize),2);
            %estimatedSize, cs
            [estimatedFlowCounts0]=PredictFlowCount(D(ii),sketchLetSize,2);
            %fill
           estimatedFlowCounts(ii) = estimatedFlowCounts0;
        end
        
        
        end
    
       %cluster + cs
    if choice ==-3
        
        %config
        %clusterMax = 4;
        
        choiceCluster = 4; %k-means
        %choiceCluster = 4; %k-means
        %cluster
        [idx,lengthPerCluster]=clusterFlow4Groups(D,choiceCluster,clusterMax);
        estimatedFlowCounts=zeros(N,1);
        %construct by flow distributions
        for i=unique(idx)'
            %index of cluster i
            ii=find(idx==i);
            %build sketch
            sketchLetSize = max(ceil((length(ii)/sum(lengthPerCluster))*sketchSize),2);
            %estimatedSize, cs
            [estimatedFlowCounts0]=PredictFlowCount(D(ii),sketchLetSize,1);
            %fill
           estimatedFlowCounts(ii) = estimatedFlowCounts0;
        end
        
        
    end
    %cluster + rps
    if choice ==-2
        
        %config
        %clusterMax = 4;
        choiceCluster = 3; %k-means
         %cluster
        [idx,lengthPerCluster]=clusterFlow4Groups(D,choiceCluster,clusterMax);
        estimatedFlowCounts=zeros(N,1);
        %construct by flow distributions
        for i=unique(idx)'
            %index of cluster i
            ii=find(idx==i);
            %build sketch
            sketchLetSize = max(ceil((length(ii)/sum(lengthPerCluster))*sketchSize),2);
            %estimatedSize
            [estimatedFlowCounts0]=PredictFlowCount(D(ii),sketchLetSize,0);
            %fill
           estimatedFlowCounts(ii) = estimatedFlowCounts0;
        end
        
        
    end
        

if choice == -1
    %separate flows
    
   percentHH = 0.1;
    counts = ceil(percentHH*N);
    [nonHH,HHs]=filterHeavyHitters(D,counts);
    
   groundtruthFlowCounts = [HHs' nonHH'];
   
   [estimatedFlowCounts0, groundtruthFlowCounts0]=PredictFlowCount(nonHH,sketchSize,0,clusterMax, pSampled,useWeight,DTrained);
   
   estimatedFlowCounts=[HHs' estimatedFlowCounts0'];
    
   return; 
end
    
    
if choice ==0
    
    %one hash
    NumNonZeros = 1;

%     U = sparse(length(D),d);
%   
%  for i=1:N
%     %select position? 
    
    signValue = 1;
%     U0 = zeros(N,1);
%      for i=1:N
%         %select position?    
%          %for indexPos= 1:NumNonZeros
% 
%             pos = ceil(rand*d);
% 
%             %signal
%             if 0
%             sigNal=rand<0.5;
%             if sigNal
%                 signValue = -1;
%             else
%                 signValue = 1;
%             end
%             end        
%             U0(i)= pos;       
%          %end    
%     end
  U0 = randi(d,[N,1]);
%sparse  
U = sparse(1:1:N,U0, signValue);

    
%     signValue = 1;
%     
%      parfor indexPos= 1:NumNonZeros
%         
%         pos = ceil(rand*d);
% 
%         %signal
%         if 0
%         sigNal=rand<0.5;
%         if sigNal
%             signValue = -1;
%         else
%             signValue = 1;
%         end
%         end
%         
%         U(i,pos)=signValue;
%     
%         
%      end
%     
%  end
 
 [ur,uc] = size(U);
 [dr,dc] = size(D);
 fprintf(1,'U:(%d %d), D:(%d %d)\n',ur,uc,dr,dc);
 XX = U'*D;

 %use the reverse
 reverse = 2;
 %svd based mp-inverse of diagonal matrix
if reverse==1
    [MPAResult]=inverseSingularMatrix(U'*U);
    inverseU = (MPAResult*U')';
end
%simply reverse the diagonal items
%worse than the mp-inverse of diagonal matrix
if reverse ==2
   [MPAResult]=inverseSingularMatrixSimple(U'*U); 
    inverseU = (MPAResult*U')';
end


%inverseU = (mldivide(U'*U,U'))';
 
  %approximated
 estimatedFlowCounts = inverseU*(XX);
  
    
end

%count sketch
if choice ==1
    
    NumNonZeros = 1;
    U = sparse(N,d); 
 
 for i=1:N
    %select position? 
    
     for indexPos= 1:NumNonZeros
        
        pos = ceil(rand*d);
        
        %signal
        sigNal=rand<0.5;
        if sigNal
            signValue = -1;
        else
            signValue = 1;
        end
        U(i,pos)=signValue;
    
        
    end
 end
 
 XX = U'*D;
 %select median
inverseU = U;
 
  %approximated
 estimatedFlowCounts  = inverseU*(XX);
end

%count min
if choice ==2
     U = sparse(length(D),d);
   
     signValue = 1;
 
 for i=1:N
    %select position? 
    
     %for indexPos= 1:NumNonZeros
        
        pos = ceil(rand*d);
   
        U(i,pos)=signValue;
    
        
    %end
 end
 
 XX = U'*D;
 %select median
inverseU = U;
 
  %approximated
  estimatedFlowCounts = inverseU*(XX);
    
end


%count sketch
if choice ==3
    
     NumOfBanks = 3; % by Elastic sketch
      
    rec=[];
    for iB=1:NumOfBanks
      
        NumNonZeros = 1;
        U = sparse(N,d); 
 
         for i=1:N
            %select position? 

             for indexPos= 1:NumNonZeros

                pos = ceil(rand*d);

                %signal
                sigNal=rand<0.5;
                if sigNal
                    signValue = -1;
                else
                    signValue = 1;
                end
                U(i,pos)=signValue;


            end
         end

         XX = U'*D;
         %select median
        inverseU = U;

          %approximated
         estimatedFlowCounts  = inverseU*(XX);
         
        rec=[rec;estimatedFlowCounts'];
    end%end
      [estimatedFlowCounts]=median(rec);
      estimatedFlowCounts=estimatedFlowCounts';
 

end

%count min
if choice ==4
    
    NumOfBanks = 3; % by Elastic sketch
    rec=[];
   
    for iB=1:NumOfBanks
    
         U = sparse(length(D),d);

         signValue = 1;

         for i=1:N
            %select position? 

             %for indexPos= 1:NumNonZeros

                pos = ceil(rand*d);

                U(i,pos)=signValue;


            %end
         end

         XX = U'*D;
         %select median
        inverseU = U;         
        estimatedFlowCounts = inverseU*(XX);

        rec=[rec;estimatedFlowCounts'];
    end%end
      [estimatedFlowCounts,idxBank]=min(rec);
      
      estimatedFlowCounts=estimatedFlowCounts';
    clear rec;
end

%fit theorem
if choice ==5
    
    NumOfBanks = 1; % by Elastic sketch
    rec=[];
   
    for iB=1:NumOfBanks
    
         U = sparse(length(D),d);

         signValue = 1;

         for i=1:N
            %select position? 

             %for indexPos= 1:NumNonZeros

                pos = ceil(rand*d);

                U(i,pos)=signValue;


            %end
         end

         XX = U'*D;
         %select median
        inverseU = U;         
        estimatedFlowCounts = inverseU*(XX);

        rec=[rec;estimatedFlowCounts'];
    end%end
      [estimatedFlowCounts,idxBank]=min(rec);
      
      estimatedFlowCounts=estimatedFlowCounts';
    clear rec;
end


fclose(ExpLog);

