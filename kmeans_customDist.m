% variation of kmeans with a custom metric "custom_metric"
function [assignment,centers] = kmeans_customDist(X,k, centers = 0, maxIter = 20000, custom_metric)
  if(centers == 0)
    centerRows = randperm(size(X)(1));
    centers = X(centerRows(1:k),:);
  endif
  disp('centers');
  disp(centers);
  #centers
  numOfRows = length(X(:,1));
  numOfFeatures = length(X(1,:));
  assignment = ones(1,numOfRows);
  i=0;
  for iter = 1:maxIter
     clusterTotals = zeros(k,numOfFeatures);
     clusterSizes = zeros(k,1);
     for rowIx = 1:numOfRows
      minDist = Inf;  ##realmax
      assignTo = 0;
      for centerIx = 1:k
        #X(rowIx,:)
        #Scenters(centerIx,:)
        #dist = sqrt(sum((X(rowIx,:)-centers(centerIx,:)).^2));
        dist = custom_metric(X(rowIx,:),centers(centerIx,:));
        if any(isnan(dist))
          disp(X(rowIx,:));
          disp(centers(centerIx,:));
          assert(1 > 1);
        endif
        if dist<minDist
          minDist = dist;
          assignTo = centerIx;
        endif
      endfor
      assignment(rowIx) = assignTo;
      clusterTotals(assignTo,:) += X(rowIx,:);
      clusterSizes(assignTo)++;
     endfor
     newCenters = zeros(k,numOfFeatures);
     for centerIx = 1:k
      newCenters(centerIx,:) = clusterTotals(centerIx,:)/clusterSizes(centerIx);
     endfor
     #dif = sqrt(sum((newCenters-centers).^2));
     ##disp('size of newCenters')
     ##disp(size(newCenters));
     ##disp('newCenters');
     ##disp(newCenters);
     ##disp('centers');
     ##disp(centers);
     ##assert(1 > 1);
     ##dif = custom_metric(newCenters,centers);
     #if dif<eps
     # break;
     #endif
     centers = newCenters;
     i++;
     ##disp(i);
     if any(isnan(centers))
       assert(2 > 2);
     endif
  endfor
  assignment= assignment';
  i
 endfunction
