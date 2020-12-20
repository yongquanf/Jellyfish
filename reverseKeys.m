function [reversedD]=reverseKeys(D,index)
%reverse key
N = max(index);
reversedD = zeros(N,1);
for i=1:length(D)
    %find the reversed index of current element; append to the last
    reversedD(index(i))=reversedD(index(i)) + D(i);
end
