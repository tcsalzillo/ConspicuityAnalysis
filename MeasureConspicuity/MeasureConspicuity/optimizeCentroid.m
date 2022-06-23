function[segCentroidX,segCentroidY] = optimizeCentroid(regInner,segCentroidX,segCentroidY,margin)
rayLength = 1:50;
%numAngles should be high enough to find a new appropriate centroid, but too
%many search angles can result in the algorithm determining that the optimal
%placement be in a narrow area that happens to span the most of regInner
%(such as the anterior parotids). 18-36 angles seems to be ideal to
%prevent this from happening consistently.
numAngles = 24; 
numTries = 3;
 
checkBoundaryX = segCentroidX-margin:segCentroidX+margin;
checkBoundaryY = segCentroidY-margin:segCentroidY+margin;
outsidePoints = find(regInner(checkBoundaryY,checkBoundaryX) == 0);
       
k = 0;
while isempty(outsidePoints) == 0
    k = k+1;
    if k>numTries
        break
    end

    rayAngles = linspace(0,(numAngles-1)*(360/numAngles),numAngles);
    rayMask = zeros(size(regInner,1),size(regInner,2)); %reset each angle
    totalRayMask = zeros(size(regInner,1),size(regInner,2)); %saved for all angles
    rayX = zeros(size(rayLength,2),size(rayAngles,2)); 
    rayY = zeros(size(rayLength,2),size(rayAngles,2));
    
    thickness = zeros(numAngles,1);
    nearestPointX = zeros(numAngles,1);
    farthestPointX = zeros(numAngles,1);
    nearestPointY = zeros(numAngles,1);
    farthestPointY = zeros(numAngles,1);
    midpointX = zeros(numAngles,1);
    midpointY = zeros(numAngles,1);

    for j = 1:size(rayAngles,2) 
        rayMask = 0*rayMask;
        rayX(:,j) = segCentroidX + round(rayLength*cosd(rayAngles(j)));
        rayY(:,j) = segCentroidY + round(rayLength*sind(rayAngles(j)));

        rayElements = (rayX(:,j)-1)*size(rayMask,1)+rayY(:,j); 
        rayMask(rayElements) = 1;
        totalRayMask(rayElements) = 1;

        centroidToShell = find((regInner(:,:) > 0 & rayMask(:,:) > 0)); 
        if isempty(centroidToShell) == 0
            xFind = floor(centroidToShell/size(regInner,1)) + 1; 
            yFind = centroidToShell - (xFind-1)*size(regInner,1); 
            mag2 = (xFind-segCentroidX).*(xFind-segCentroidX) + (yFind-segCentroidY).*(yFind-segCentroidY);

            [~, minMag] = min(mag2);
            [~, maxMag] = max(mag2);

            nearestPointX(j) = xFind(minMag);
            farthestPointX(j) = xFind(maxMag);
            nearestPointY(j) = yFind(minMag);
            farthestPointY(j) = yFind(maxMag);

            thickness(j) = sqrt(mag2(maxMag)) - sqrt(mag2(minMag));
            midpointX(j) = round((farthestPointX(j)+nearestPointX(j))/2);
            midpointY(j) = round((farthestPointY(j)+nearestPointY(j))/2);

            tempStart = find(rayElements == centroidToShell(minMag));
            tempStop = find(rayElements == centroidToShell(maxMag));
            if any(regInner(rayElements(tempStart:tempStop)) == 0)
                thickness(j) = -999;
            end

            checkNewBoundaryX = midpointX(j)-margin:midpointX(j)+margin;
            checkNewBoundaryY = midpointY(j)-margin:midpointY(j)+margin;
            newOutsidePoints = find(regInner(checkNewBoundaryY,checkNewBoundaryX) == 0);
            if (regInner(segCentroidY,segCentroidX)>0) && (length(newOutsidePoints) >= length(outsidePoints))
                thickness(j) = -888;
            end
        end
    end   

    [maxThickness, maxThicknessIndex] = max(thickness(:));
    if maxThickness < 1
        break
    end
    
    segCentroidX = midpointX(maxThicknessIndex);
    segCentroidY = midpointY(maxThicknessIndex);
       
    checkBoundaryX = segCentroidX-margin:segCentroidX+margin;
    checkBoundaryY = segCentroidY-margin:segCentroidY+margin;
    outsidePoints = find(regInner(checkBoundaryY,checkBoundaryX) == 0);

%     bestRayMask = 0*rayMask;
%     bestRayX(:) = segCentroidX + round(rayLength*cosd(rayAngles(maxThicknessIndex)));
%     bestRayY(:) = segCentroidY + round(rayLength*sind(rayAngles(maxThicknessIndex)));
%     bestRayElements = (bestRayX(:)-1)*size(bestRayMask,1)+bestRayY(:);
% 
%     bestRayMask(regInner>0) = 75;
%     bestRayMask(bestRayElements) = 100;
%     bestRayMask(segCentroidY-2:segCentroidY+2,segCentroidX) = 150;
%     bestRayMask(segCentroidY,segCentroidX-2:segCentroidX+2) = 150;
% 
%     bestRayMask(segCentroidY+nearestPointY(maxThicknessIndex),segCentroidX+nearestPointX(maxThicknessIndex)) = 300;
%     bestRayMask(segCentroidY+farthestPointY(maxThicknessIndex),segCentroidX+farthestPointX(maxThicknessIndex)) = 300;
% 
%     bestRayMask(segCentroidY-2:segCentroidY+2,segCentroidX) = 300;
%     bestRayMask(segCentroidY,segCentroidX-2:segCentroidX+2) = 300;
%     bestRayMask(nearestPointY(maxThicknessIndex),nearestPointX(maxThicknessIndex)) = 300;               
%     bestRayMask(farthestPointY(maxThicknessIndex),farthestPointX(maxThicknessIndex)) = 300;
% 
%     tempfig2 = figure;
%     imagesc(bestRayMask)
%     colormap gray
%     uiwait(tempfig2)
end



