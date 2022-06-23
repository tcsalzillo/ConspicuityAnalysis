function[conspicuityData,flag] = measureConspicuity(dicomStack,segVolume,segShell,segVolumeDicoms,segShellDicoms)
%This function calculates conspicuity based on paper DOI:10.1097/00004424-197411000-00009 
%Essentially, the script loops through each slice and determines how many
%unique non-connected regions in the structure are present. For each 
%region, the centroid of the innerSegVolume is calculated (or the average
%position if there are multiple subregions in the particular region). 
%
%This centroid serves as the source of several rays (set by numAngles)
%which expand radially outwards (set by rayLength) and whose direction 
%spans from zero to 360 degrees in increments of 360/numAngles. In order
%To determine each unique set of inner, struct, and outer points of the
%segmentation, used in the conspicuity measurement, the script cacluates
%the intersection points of each ray and the inner, struct, and outer 
%segShells. 
%
%Any ray that does not intersect with each segShell is removed.
%Checks are in place to make sure that there is only one innerSegShell
%pixel for each structSegShell and outerSegShell intersect point. In 
%In other words, pixels within each SegShell are used only once 
%in the total set of intersection points. This ensures 1:1:1 matching in
%the conspicuity measurement. 

close all
conspicuityData = {'NA','NA','NA','NA';'NA','NA','NA','NA'};
flag = 0;
displayFig = 1; %Set to 0 if not wanting to view figures

fprintf('Measuring Conspicuity\n')

numChecks = 10;
rayLength = [0:200]; %Set length to include all of segmentation
numRotations = 1;
numAngles = 360; %Set to number of angles you wish to interrogate
%You really only need more than one rotation when numAngles is small. The
%point of using multiple rotations is to capture any inner-struct-outer
%shell triplets that were missed when spanning the numAngles. Furthermore,
%we only want one outer point for each inner point, so as long as each of
%the innerShell points are traversed, the total number of accepted rays won't
%change. 

sliceLabel = [];
regionLabel = [];
conspicuity = [];
counter = 0; %increases during each loop to add to conspicuity vector

%Determine which slices have segmentations present
segSlice = zeros(size(segShell,3),1); %set of slices with a segmentation present
for j=1:size(segShell,3)
    if any(reshape(segShell(:,:,j,3).',1,[])) % Determines if slice contains
        %any values greater than zero (if segmentation is present)
        segSlice(j) = j;
    end
end
segSlice(segSlice == 0) = [];
numSegSlice = size(segSlice,1);

if isempty(segSlice)
    %fprintf('\nSegmentation is empty. Moving to next Observer.\n')
    flag = 1;
    return
end

if displayFig == 1
    %Determine the appropriate viewing window for image display
    tempfig = figure;
    segMIP = sum(segVolume(:,:,:,1),3); %sums each slice to see structure extent
    segMIP(find(segMIP)) = 1; %creates maximum intensity projection (MIP)
    xDimMIP = sum(segMIP,1);
    yDimMIP = sum(segMIP,2);
    %below sets viewing window to encompass MIP (and thus the segmentation on
    %each slice) with a 20 pixel border
    imageRangeX = find(xDimMIP,1,'first')-20:find(xDimMIP,1,'last')+20;
    imageRangeY = find(yDimMIP,1,'first')-20:find(yDimMIP,1,'last')+20;
    %Check to make sure imageRange is within image, otherwise adjust
    imageRangeX(imageRangeX < 1) = [];
    imageRangeX(imageRangeX > size(dicomStack,2)) = [];
    imageRangeY(imageRangeY < 1) = [];
    imageRangeY(imageRangeY > size(dicomStack,1)) = [];
    imagesc(segMIP(imageRangeY,imageRangeX)) 
    colormap gray
    movegui('northwest')
    uiwait(tempfig,2) %pauses program to check viewing window 
                    %comment out to skip the pause
%     fprintf('Check image range \n') %place conditional breakpoint here

    tempfig1 = figure;
    tempfig2 = figure;
end

fprintf('Final Slice: %3.0f ; Current Slice: %3.0f',segSlice(numSegSlice), 0)
for g = 1:numSegSlice
    slice = segSlice(g);
    fprintf('\b\b\b%3.0f',slice)
    
%   Find unique regions in slice
    innerBW = segVolume(:,:,slice,3) > 0;
    innerSegStats = regionprops(bwconncomp(innerBW),'all'); 
    structBW = segVolume(:,:,slice,2) > 0;
    structSegStats = regionprops(bwconncomp(structBW),'all');
    outerBW = segVolume(:,:,slice,1) > 0;
    outerSegStats = regionprops(bwconncomp(outerBW),'all'); 

    numRegion = min([size(innerSegStats,1),size(structSegStats,1),size(outerSegStats,1)]); %determines number of unique segmentations on slice

    for h = 1:numRegion %allows you to analyze each unique segmentation on slice 
    %for h = 2:2
        rotConspicuity = zeros(numRotations,1);
        deletePoint = zeros(numRotations,1);       
        trackedOuterShellPoints = zeros(numAngles,numRotations); %to find duplicate outerShellPoints

        regInner = zeros(size(dicomStack,1),size(dicomStack,2));
        regStruct = zeros(size(dicomStack,1),size(dicomStack,2));
        regOuter = zeros(size(dicomStack,1),size(dicomStack,2));
        regInner(innerSegStats(h).PixelIdxList) = 1; %Extracts seg pixels from unique segmentation
        regStruct(structSegStats(h).PixelIdxList) = 1;
        regOuter(outerSegStats(h).PixelIdxList) = 1;

        segCentroid = round(innerSegStats(h).Centroid); %rounds to nearest pixel
        segCentroidX = segCentroid(1);
        segCentroidY = segCentroid(2);

        %Determine how close centroid is to segmentation boundary by
        %calculating number of surrounding pixels that equal zero
        margin = 3;
        checkBoundaryX = segCentroidX-margin:segCentroidX+margin;
        checkBoundaryY = segCentroidY-margin:segCentroidY+margin;
        outsidePoints = find(regInner(checkBoundaryY,checkBoundaryX) == 0);
        
        %Moves ray source from centroid if too close to edge.
        if isempty(outsidePoints) == 0
            centroidMask = 0*regInner;
            centroidMask(edge(regInner)>0) = 200;
            centroidMask(segCentroidY-2:segCentroidY+2,segCentroidX) = 100;
            centroidMask(segCentroidY,segCentroidX-2:segCentroidX+2) = 100;

            %Will move the source of rays to a more optimal location (away
            %from segmentation boundaries) if enabled
            [segCentroidX,segCentroidY] = optimizeCentroid(regInner,segCentroidX,segCentroidY,margin);
%                 fprintf('\n centroid updated on slice %d',slice);

            centroidMask(segCentroidY-2:segCentroidY+2,segCentroidX) = 300;
            centroidMask(segCentroidY,segCentroidX-2:segCentroidX+2) = 300;

            if displayFig == 1
                figure(tempfig1);
                imagesc(centroidMask(imageRangeY,imageRangeX))
                colormap gray
                movegui('northeast')
    %             uiwait(tempfig1,1.5)
                drawnow
            end
        end

        for p = 1:numRotations
            %rayAngles = linspace((0+((-1)^p)*(p-1)*(360/numRotations)),((-1)^p)*(((numAngles-1)*(360/numAngles)+(p-1)*(360/numRotations))),numAngles); 
            %alters starting angle to span all ray intersect angles and 
            %Doesn't count zero and 360 ray twice
            rayAngles = linspace(0,(numAngles-1)*(360/numAngles),numAngles); %Don't count 
            %Zero and 360 ray twice            
            if mod(p,2) == 0
                rayAngles = linspace(360,0+(numAngles-1)*(1/numAngles),numAngles); %swaps direction the angle vector is 
                %spanned to capture all intersect points
            end

            rayMask = zeros(size(dicomStack,1),size(dicomStack,2)); %reset each angle
            totalRayMask = zeros(size(dicomStack,1),size(dicomStack,2)); %saved for all angles
            rayX = zeros(size(rayLength,2),size(rayAngles,2)); 
            rayY = zeros(size(rayLength,2),size(rayAngles,2));
            innerShellPoints = zeros(size(rayAngles,2),1); %intersection points
            structShellPoints = zeros(size(rayAngles,2),1);
            outerShellPoints = zeros(size(rayAngles,2),1);
            removePoints = zeros(size(rayAngles,2),1); %tag to remove intersection points          

            %below if statement handles case where there is 1 outer seg shell,
            %but multiple inner shells (likely due to the structure being small
            %in one dimension, so the edge function splits it into multiple
            %inner shells). Thus the mean location of each inner shell centroid
            %is determined and all inner shell regions are considered in the
            %ray intersect algorithm (instead of just the first region)
            if ((size(outerSegStats,1) == 1) && (size(innerSegStats,1) > 1))
                totalCentroidX = 0;
                totalCentroidY = 0;
                for t = 1:size(innerSegStats,1)
                    totalCentroidX = totalCentroidX + innerSegStats(t).Centroid(1);
                    totalCentroidY = totalCentroidY + innerSegStats(t).Centroid(2);
                    regInner(innerSegStats(t).PixelIdxList) = 1;
                end

                for t = 1:size(structSegStats,1)
                    regStruct(structSegStats(t).PixelIdxList) = 1;
                end        

                segCentroidX = round(totalCentroidX/size(innerSegStats,1));
                segCentroidY = round(totalCentroidY/size(innerSegStats,1));
             end

    %         figure
    %         imagesc(regInner(imageRangeY,imageRangeX))
    %         colormap gray

            for j = 1:size(rayAngles,2) %This loop determines the intersection points of the 
                %rays and various segmentation shells
                rayMask = 0*rayMask;
                
                rayX(:,j) = segCentroidX + round(rayLength*cosd(rayAngles(j)));
                rayY(:,j) = segCentroidY + round(rayLength*sind(rayAngles(j)));      
                
                tempRayX = rayX(:,j);
                tempRayY = rayY(:,j);
                
                tempRayX(rayX(:,j)>size(dicomStack,2)) = [];
                tempRayY(rayX(:,j)>size(dicomStack,2)) = [];
                tempRayX(rayY(:,j)>size(dicomStack,1)) = [];
                tempRayY(rayY(:,j)>size(dicomStack,1)) = [];
                
                rayElements = (tempRayX(:)-1)*size(rayMask,1)+tempRayY(:); %Converts X and Y ray             
                %vectors to the element number in matrix (matrix is column-major)and
                rayMask(rayElements) = 1;
                totalRayMask(rayElements) = 1;

                %Below is a calculation to find the outermost intersection
                %point between ray and seg shell. This is only necessary when
                %shell is more than 1 pixel thick (such as at top and bottom
                %slice of segmentation). The problem of using 'last' (as in the
                %commented command) is that it returns the numerically largest
                %index of intersection. For rays whose indeces decrease as they
                %move outward (those moving leftwards and downwards), this
                %results in the innermost intersection point being selected,
                %not the outermost.
                %tempInnerShell = find((regInner(:,:) > 0 & rayMask(:,:) > 0),1,'last');
                tempInnerShell = find((regInner(:,:) > 0 & rayMask(:,:) > 0)); %these are all intersetction points
                xFind = floor(tempInnerShell/size(dicomStack,1)) + 1; %calculates x coordinate of intersection point
                yFind = tempInnerShell - (xFind-1)*size(dicomStack,1); %calculates y coordinate of intersection point
                mag2 = (xFind-segCentroidX).*(xFind-segCentroidX) + (yFind-segCentroidY).*(yFind-segCentroidY);
                %Calculates displacement between centroid and intersection
                %point
                tempInnerShell(mag2 ~= max(mag2)) = []; %deletes all intersection points less than max displacement

                if isempty(tempInnerShell) %flag if ray does not intersect shell (usually between diagonally-connected points)
                    removePoints(j) = 1; 
                    innerShellPoints(j) = -999; %ensures point is not included
    %                 fprintf('Ray %d does not intersect with inner segShell \n',rayAngles(j))
                else
                    innerShellPoints(j) = tempInnerShell;
                end

                tempStructShell = find((regStruct > 0 & rayMask(:,:) > 0));
                xFind = floor(tempStructShell/size(dicomStack,1)) + 1;
                yFind = tempStructShell - (xFind-1)*size(dicomStack,1);
                mag2 = (xFind-segCentroidX).*(xFind-segCentroidX) + (yFind-segCentroidY).*(yFind-segCentroidY);
                tempStructShell(mag2 ~= max(mag2)) = [];

                if isempty(tempStructShell)
                    removePoints(j) = 1;
                    structShellPoints(j) = -999;
                    innerShellPoints(j) = -999; %ensures this point is not considered as duplicate
                      % if not included, a point not flagged in inner but
                      % flagged in struct or outer will be considered for
                      % duplication
    %                 fprintf('Ray %d does not intersect with struct segShell \n',rayAngles(j))
                else
                    structShellPoints(j) = tempStructShell;
                end

                tempOuterShell = find((regOuter > 0 & rayMask(:,:) > 0));
                xFind = floor(tempOuterShell/size(dicomStack,1)) + 1;
                yFind = tempOuterShell - (xFind-1)*size(dicomStack,1);
                mag2 = (xFind-segCentroidX).*(xFind-segCentroidX) + (yFind-segCentroidY).*(yFind-segCentroidY);
                tempOuterShell(mag2 ~= max(mag2)) = [];

                if isempty(tempOuterShell)
                    removePoints(j) = 1;
                    innerShellPoints(j) = -999;
                    outerShellPoints(j) = -999; %ensures this point is not considered as duplicate
    %                 fprintf('Ray %d does not intersect with outer segshell \n',rayAngles(j))
                else
                    outerShellPoints(j) = tempOuterShell;
                    trackedOuterShellPoints(j,p) = tempOuterShell;
                end
                
                
%                 tempMask = regInner+regStruct+regOuter;
%                 tempMask(tempInnerShell) = 5;
%                 tempMask(tempStructShell) = 5;
%                 tempMask(tempOuterShell) = 5;
%                 figure
%                 imagesc(tempMask)

                %Now that the rays that intersect *each shell* have been
                %identified, identify the rays that produce duplicate
                %intersection points with the inner shell and remove. In other
                %words, each iteration of ray angles determines if it produces
                %intersection points that are identical to any of the prior
                %iterations
                if (j > 1) && (isempty(tempInnerShell) == 0) && (any(find(innerShellPoints(1:j-1) == tempInnerShell)))  
                    removePoints(j) = 1;
    %                 fmt = ['Ray %d intersects at %d which already exists at rays' ... 
    %                     repmat(' %1.0f ',1,numel(rayAngles(find(innerShellPoints(1:j-1) == tempInnerShell)))) '\n'];
    %                 fprintf(fmt,rayAngles(j),tempInnerShell,rayAngles(find(innerShellPoints(1:j-1) == tempInnerShell)))
                    %sprintf('%d%s%d',j,' ; ',tempInnerShell)
                    %1. Repeat pixels in inner should propogate to struct and outer
                    %since their diameters are larger
                    %2. Implementing this on each structure for large number of angles
                    %makes it nearly impossible for each ray to hit a unique pixel for
                    %each structure- this is why only the inner struct is
                    %considered
                    innerShellPoints(j) = -888; %ensures points are not included
                    structShellPoints(j) = -888;
                    outerShellPoints(j) = -888;
                end

                %below checks if any of the outershell points in rotation p
                %match the outershell points in p-1. It can be assummed
                %that two rays originating at the centroid intersecting
                %with the same outershell point are duplicate and should
                %not be muliplty counted (otherwise they are weighted more
                %than those rays captured in only one rotation scheme)
                if((p>1) && (isempty(tempOuterShell) == 0) && (any(find(trackedOuterShellPoints(:,1:p-1) == tempOuterShell))))
                    removePoints(j) = 1;
                    innerShellPoints(j) = -777; %ensures points are not included
                    structShellPoints(j) = -777;
                    outerShellPoints(j) = -777;
                end

                %Below removes points if the shells "overlap" (if the
                %intersection point in one shell is the same location as
                %the intersection point in another shell). These points
                %should be nested, not overlapping, but the grow/shrink
                %command is not perfect
                if ((isequal(tempInnerShell,tempStructShell)) || (isequal(tempStructShell,tempOuterShell)) || (isequal(tempInnerShell,tempOuterShell)))
                    removePoints(j) = 1;
                    
%                     sprintf('%d%s%d%s%d',tempInnerShell,' ',tempStructShell,' ',tempOuterShell)
                    
                    innerShellPoints(j) = -555;
                    structShellPoints(j) = -555;
                    outerShellPoints(j) = -555;
                end
            end

            %figure of original seg shells
            if displayFig == 1
                segShellFigure = zeros(size(dicomStack,1),size(dicomStack,2));
                segShellFigure(segShell(:,:,slice,1)>0) = 150;
                segShellFigure(segShell(:,:,slice,2)>0) = 100;
                segShellFigure(segShell(:,:,slice,3)>0) = 50;
                segShellFigure(segCentroidY-2:segCentroidY+2,segCentroidX) = 150;
                segShellFigure(segCentroidY,segCentroidX-2:segCentroidX+2) = 150;

                totalRayMask = zeros(size(dicomStack,1),size(dicomStack,2));
                totalRayMask(segShell(:,:,slice,1)>0) = 75;
                totalRayMask(segShell(:,:,slice,2)>0) = 75;
                totalRayMask(segShell(:,:,slice,3)>0) = 75;
                for j=1:size(rayAngles,2)
                    tempRayX = rayX(:,j);
                    tempRayY = rayY(:,j);

                    tempRayX(rayX(:,j)>size(dicomStack,2)) = [];
                    tempRayY(rayX(:,j)>size(dicomStack,2)) = [];
                    tempRayX(rayY(:,j)>size(dicomStack,1)) = [];
                    tempRayY(rayY(:,j)>size(dicomStack,1)) = [];

                    totalRayMask((tempRayX(:)-1)*size(rayMask,1)+tempRayY(:)) = 100*(1.5-0.5*(-1)^j); 
                   %Visualizes rays with alternating brightness for clarity
                end
        %         totalRayMask(outerShellPoints) = totalRayMask(outerShellPoints) + 50;
        %         totalRayMask(structShellPoints) = totalRayMask(structShellPoints) + 50;
        %         totalRayMask(innerShellPoints) = totalRayMask(innerShellPoints) + 50;
                totalRayMask(segCentroidY-2:segCentroidY+2,segCentroidX) = 250;
                totalRayMask(segCentroidY,segCentroidX-2:segCentroidX+2) = 250;

        %         figure 
        %         imagesc(totalRayMask(imageRangeY,imageRangeX))
        %         colormap gray
            end

            % sprintf(['Before edit, outerShellPoints =' repmat(' %1.0f',1,numel(outerShellPoints))],outerShellPoints)
            % sprintf(['Before edit, structShellPoints =' repmat(' %1.0f',1,numel(structShellPoints))],structShellPoints)
            % sprintf(['Before edit, innerShellPoints =' repmat(' %1.0f',1,numel(innerShellPoints))],innerShellPoints)

            %Check Ray Angle vectors that were valid
            %tempRayAngles = rayAngles;
            %tempRayAngles(removePoints == 1) = -999;
            %rayAngleDocument(:,size(rayAngleDocument,2)+1) = cat(1,slice,h,p,tempRayAngles(:));

            %Remove data whose rays were flagged
            rayAngles(removePoints == 1) = [];
            rayX(:,removePoints == 1) = [];
            rayY(:,removePoints == 1) = [];

            outerShellPoints(removePoints == 1) = [];
            structShellPoints(removePoints == 1) = [];
            innerShellPoints(removePoints == 1) = [];

            finalNumAngles = size(rayAngles,2);
    %         fprintf('Final num angles = %d \n', finalNumAngles)

            % sprintf(['After edit, outerShellPoints =' repmat(' %1.0f',1,numel(outerShellPoints))],outerShellPoints)
            % sprintf(['After edit, structShellPoints =' repmat(' %1.0f',1,numel(structShellPoints))],structShellPoints)
            % sprintf(['After edit, innerShellPoints =' repmat(' %1.0f',1,numel(innerShellPoints))],innerShellPoints)

            if displayFig == 1
                finalRayMask = zeros(size(dicomStack,1),size(dicomStack,2));
                finalRayMask(segShell(:,:,slice,1)>0) = 75;
                finalRayMask(segShell(:,:,slice,2)>0) = 75;
                finalRayMask(segShell(:,:,slice,3)>0) = 75;
                for j=1:finalNumAngles
                    tempRayX = rayX(:,j);
                    tempRayY = rayY(:,j);

                    tempRayX(rayX(:,j)>size(dicomStack,2)) = [];
                    tempRayY(rayX(:,j)>size(dicomStack,2)) = [];
                    tempRayX(rayY(:,j)>size(dicomStack,1)) = [];
                    tempRayY(rayY(:,j)>size(dicomStack,1)) = [];                    
                    
                    finalRayMask((tempRayX(:)-1)*size(rayMask,1)+tempRayY(:)) = 100*(1.5-0.5*(-1)^j); 
                end
                finalRayMask(outerShellPoints) = finalRayMask(outerShellPoints) + 50;
                finalRayMask(structShellPoints) = finalRayMask(structShellPoints) + 50;
                finalRayMask(innerShellPoints) = finalRayMask(innerShellPoints) + 50;
                finalRayMask(segCentroidY-2:segCentroidY+2,segCentroidX) = 250;
                finalRayMask(segCentroidY,segCentroidX-2:segCentroidX+2) = 250;
            end

            if isempty(rayAngles) == 1
                deletePoint(p) = 1;
            else   
                tempOuter = segVolumeDicoms(:,:,slice,1);
                tempStruct = segVolumeDicoms(:,:,slice,2);
                tempInner = segVolumeDicoms(:,:,slice,3);
                
                if displayFig == 1
                    %Create label map that includes the seg shell values that intersect with
                    %the acceptable rays (those that intersect with each of the shells)
                    segRayIntercept = zeros(size(dicomStack,1),size(dicomStack,2));
                    segRayIntercept(outerShellPoints) = 150;
                    segRayIntercept(structShellPoints) = 100;
                    segRayIntercept(innerShellPoints) = 50;
                    
                    figSegRayIntDicom = zeros(size(dicomStack,1),size(dicomStack,2));
                    figSegRayIntDicom(outerShellPoints) = tempOuter(outerShellPoints);
                    figSegRayIntDicom(structShellPoints) = tempStruct(structShellPoints);
                    figSegRayIntDicom(innerShellPoints) = tempInner(innerShellPoints);
   
%                     figure
%                     imagesc(segRayIntercept(imageRangeY,imageRangeX))
%                     colormap gray
%                     figure
%                     imagesc(figSegRayIntDicom(imageRangeY,imageRangeX))
%                     colormap gray
                end

                %Below creates the vector that inputs the intersection points in an
                %ordered fashion (by polar angle)
                outerDicomVec = zeros(finalNumAngles,1);
                structDicomVec = zeros(finalNumAngles,1);
                innerDicomVec = zeros(finalNumAngles,1);

                outerDicomVec = tempOuter(outerShellPoints);
                structDicomVec = tempStruct(structShellPoints);
                innerDicomVec = tempInner(innerShellPoints);

                %"Circularizes" the matrix so that j-1 and j+1 operations can be
                %performed at endpoints
                catOuterDicomVec = cat(1,outerDicomVec(end),outerDicomVec);
                catOuterDicomVec = cat(1,catOuterDicomVec,outerDicomVec(1));
                catStructDicomVec = cat(1,structDicomVec(end),structDicomVec);
                catStructDicomVec = cat(1,catStructDicomVec,structDicomVec(1));
                catInnerDicomVec = cat(1,innerDicomVec(end),innerDicomVec);
                catInnerDicomVec = cat(1,catInnerDicomVec,innerDicomVec(1));

                contrast = zeros(finalNumAngles,1);
                laplacian = zeros(finalNumAngles,1);

                %From cited papers
                for j = 2:finalNumAngles+1
                    contrast(j-1) = abs(outerDicomVec(j-1) - innerDicomVec(j-1));

                    laplacian(j-1) = (1/4)*(abs(catInnerDicomVec(j+1) + catInnerDicomVec(j-1)-...
                        2*catInnerDicomVec(j)) + abs(catOuterDicomVec(j+1) +...
                        catOuterDicomVec(j-1) - 2*catOuterDicomVec(j)));
                end

                meanContrast = mean(contrast);
                meanLaplacian = mean(laplacian);

                if meanLaplacian > 0
        %             counter = counter + 1;
        %             conspicuity(counter) = meanContrast/meanLaplacian;
        %             sliceLabel(counter) = slice;
        %             regionLabel(counter) = h;

                    rotConspicuity(p) = meanContrast/meanLaplacian;
                else
                    %sprintf('meanLaplacian = 0') %prevents divide by zero
                    deletePoint(p) = 1;
                end
            end

            %Comment out the figures below if you want to skip QA
            if ((displayFig == 1))% && (ismember(slice,checkPoints) || numRegion > 1))
                %pause(1.5) %only keep if you have >1 rotation
%                 close all
                %put breakpoint here to pause (p==1)
                figure(tempfig2) 
                imagesc(finalRayMask(imageRangeY,imageRangeX))
                colormap gray
                movegui('southeast')
                drawnow 
            end %place conditional breakpoint here (ismember(slice,checkPoints))          
        end

        rotConspicuity(deletePoint==1) = [];
        if isempty(rotConspicuity) == 0
            counter = counter + 1;
            conspicuity(counter) = mean(rotConspicuity(:));
            sliceLabel(counter) = slice;
            regionLabel(counter) = h;
        end
    end %Place conditional breakpoint here to QA each step
end    

if isempty(conspicuity)
    flag = 2;
%     fprintf('No conspicuity data measured')
    return
else
    fprintf('\nNumber of conspicuity measurements: %d \n', size(conspicuity,2))
    fprintf('Mean of total conspicuity is: %d \n',mean(conspicuity))
    fprintf('Min conspicuity: %d ; Max conspicuity: %d \n',min(conspicuity),max(conspicuity))
    
%Keep inner 90th percentile to remove noise
    tempConspicuity = sort(conspicuity);
    tempConspicuity(counter-round(0.05*counter)+1:counter) = [];
    tempConspicuity(1:round(0.05*counter)) = [];
    overallConspicuity = mean(tempConspicuity);

    t1 = {'Min Conspicuity',num2str(min(conspicuity),'%.4f'),'Max Conspicuity',num2str(max(conspicuity),'%.4f');...
        'Mean Conspicuity',num2str(mean(conspicuity),'%.4f'),'Mean Conspicuity (inner 90%)',num2str(overallConspicuity,'%.4f');...
        'Number','Slice','Region','Conspicuity'};

    t2 = [1:counter];
    t2 = cat(1,t2,sliceLabel);
    t2 = cat(1,t2,regionLabel);
    t2 = cat(1,t2,conspicuity);
    t2 = num2cell(t2');

    conspicuityData = [t1;t2];
end

if displayFig == 1
    figure
    imagesc(segShellFigure(imageRangeY,imageRangeX))
    colormap gray
    movegui('southwest')

    figure
    imagesc(segRayIntercept(imageRangeY,imageRangeX))
    %imagesc(segRayIntercept(:,:))
    colormap gray
    movegui('southwest')

    figure
    imagesc(figSegRayIntDicom(imageRangeY,imageRangeX))
    colormap gray
    movegui('southwest')

    figure 
    imagesc(totalRayMask(imageRangeY,imageRangeX))
    colormap gray
    movegui('southwest')
    
    figure 
    imagesc(finalRayMask(imageRangeY,imageRangeX))
    colormap gray
    movegui('southwest')

    figure
    histogram(conspicuity,0:0.1:8)
    %histogram(conspicuity,1:0.1:6)
    movegui('southwest')
    
    figure
    histogram(tempConspicuity,0:0.1:8)
    movegui('southwest')

    %fprintf('Final num angles = %d \n', finalNumAngles)
    pause(5);
end