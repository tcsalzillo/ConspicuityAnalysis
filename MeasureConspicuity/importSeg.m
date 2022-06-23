function[segVolume,segShell,segVolumeDicoms,segShellDicoms,flag] = importSeg(dirName,structure,observerName,dicomStack,pixShrink,pixGrow)
%This function imports individual structure segmentations located in
%\StructureX subdirectories. Segmentations should be in .nii.gz format (or
%modify the code accordingly) with only 1 segmentation present per file.
%The name of the structure segmentation should be contained in structure
%variable. This code appends the observer number provided by caseNum.
%Furthermore the last portion of the structure segmentation name is
%provided by structLabels below. Thus in this example, there should be 5
%segmentations per observer named as the following:
%StructureX_ObY - grow  2.nii.gz
%StructureX_ObY - grow  1.nii.gz
%StructureX_ObY.nii.gz
%StructureX_ObY - shrink  1.nii.gz
%StructureX_ObY - shrink  2.nii.gz
%These structLabels are used for conspicuity measurements (which is why
%there are grown and shrunk copies of each segmentation. Modify if using
%for alternate applications or have different naming conventions.

%After importing the structure segmentations, they are exported as binary
%matrices (segVolume). Then for conspicuity applications, segVolumes are
%used to create outer shells of the segmentation (segShell) which are
%exported as binary matrices. The original DICOM values are also copied 
%into segVolume and segShell locations and exported as matrices 
%(segVolumeDicoms and segShellDicoms)

%PixShrink refers to number of pixels to shrink the structure seg to produce
%the innerSegShell. By default, this is set to 2, but for small structures
%such as lymph nodes, you may only want to shrink by 1 pixel so that there
%are sufficient pixels in the innerSegShell. PixGrow can also be adjusted
%to 1 or 2, depending on how far away from the segmentation you wish to
%measure conspicuity. By default, this is set to 2 to allow for 1 mm

%May need to adjust the labels below if using a different naming convention
%structLabels = {' - grow  2.nii.gz';' - grow  1.nii.gz';'.nii.gz';' - shrink  1.nii.gz';' - shrink  2.nii.gz'};
if ((pixGrow == 1) && (pixShrink == 1))
    structLabels = {' - grow  1.nii.gz';'.nii.gz';' - shrink  1.nii.gz'};
elseif ((pixGrow == 1) && (pixShrink == 2))
    structLabels = {' - grow  1.nii.gz';'.nii.gz';' - shrink  2.nii.gz';};
elseif ((pixGrow == 2) && (pixShrink == 1))
    structLabels = {' - grow  2.nii.gz';'.nii.gz';' - shrink  1.nii.gz';};
elseif ((pixGrow == 2) && (pixShrink == 2))
    structLabels = {' - grow  2.nii.gz';'.nii.gz';' - shrink  2.nii.gz';};
else
    fprintf('Define appropriate pixel grow/shrink values')
    flag = 1;
    return
end
possibleStructLabels = {' - grow  2.nii.gz';' - grow  1.nii.gz';'.nii.gz'; ...
    ' - shrink  1.nii.gz';' - shrink  2.nii.gz'};

numLabels = size(structLabels,1);

segVolume = [];
segShell = [];
segShellDicoms = [];
segVolumeDicoms = [];
flag = 0;
fprintf('Importing Segmentations\n')

%If reading in nrrd files, uncomment below- might need to adjust subsequent
%processing
%segHeader = nhdr_nrrd_read(sprintf('%s%s', dirName, '\slicerSeg\tumor_Mod.nrrd'), 1);
%make sure nrrd_read_write_rensonnet folder is in path
% structNamePartial = sprintf('%s%s%s%s%s%d%s%s',dirName,'\',structure,'\','case_',mrn,'_',observerName);
searchPath = sprintf('%s%s%s%s',dirName,'\',structure,'\');
foundItems = dir([searchPath,'*',observerName,'*']);
if isempty(foundItems)
    flag = 1;
    return
end
segName = foundItems(1).name;

%remove structLabels from name
k=0;
clear('tempExtPos')
for j=1:size(possibleStructLabels,1)
    if isempty(strfind(segName,possibleStructLabels{j}))==0
        k = k+1;
        tempExtPos(k) = strfind(segName,possibleStructLabels{j});
    end
end
extPos = min(tempExtPos);
segName(extPos:end) = [];

structNamePartial = sprintf('%s%s',searchPath,segName);

%Check that all of the segmentation files are present
for j=1:numLabels
    if ~any(size(dir([structNamePartial char(structLabels(j,1)) ]),1))
%         fprintf('Missing sub-structure: %s (check filename)\n',char(structLabels(j,1)))
        flag = 1;
        return
    end
end

%Read in first segmentation file and initiate segVolume variable
tempSeg = niftiread(sprintf('%s%s', structNamePartial,char(structLabels(1,1))));
segVolume = zeros(size(tempSeg,1),size(tempSeg,2),size(tempSeg,3),numLabels);
segVolume(:,:,:,1) = tempSeg;
%Read in remaining segmentation files
for j=2:numLabels
    segVolume(:,:,:,j) = niftiread(sprintf('%s%s', structNamePartial,char(structLabels(j,1))));
end

%For some reason, segmentation coordinates are transposed and 
%slice order may be flipped
segVolume = permute(segVolume,[2,1,3,4]);
%segVolume = flip(segVolume,3);

%Create segmentation shells according to the order of structLabels
%Edge function does return some pixels outside of segVolume, rather than
%constraining returned pixels to segVolume area. Thus, it will only be used
%to visualize edges, while calculations use exclusively segVolume
segShell = 0*segVolume;
for j=1:size(segVolume,3)
    segShell(:,:,j,1) = edge(segVolume(:,:,j,1)); %outerShell
    segShell(:,:,j,2) = edge(segVolume(:,:,j,2)); %structShell
    segShell(:,:,j,3) = edge(segVolume(:,:,j,3)); %innerShell
end

%Below copies the original DICOM values into segVolumeDicoms and 
%segShellDicoms. Then the pixels in which the segVolume and segShell
%are less than 1 (outside segmentation) are set to zero in segVolumeDicoms
%and segShellDicoms
segVolumeDicoms = 0*segVolume;
segShellDicoms = 0*segShell;
for j = 1:3 
    tempVolumeOutput = dicomStack;
    tempVolumeOutput(segVolume(:,:,:,j)<1) = 0;
    segVolumeDicoms(:,:,:,j) = tempVolumeOutput;
end

for j = 1:3
    tempShellOutput = dicomStack;
    tempShellOutput(segShell(:,:,:,j)<1) = 0;
    segShellDicoms(:,:,:,j) = tempShellOutput;
end


% Uncomment section below to QA
% fprintf('QA starting\n')
% %First create MIP and constrain viewing window to structure
% segMIP = sum(segVolume(:,:,:,1),3); %sums each slice to see structure extent
% segMIP(find(segMIP)) = 1; %creates maximum intensity projection (MIP)
% xDimMIP = sum(segMIP,1);
% yDimMIP = sum(segMIP,2);
% %below sets viewing window to encompass MIP (and thus the segmentation on
% %each slice) with a specified pixel border
% imageRangeX = find(xDimMIP,1,'first')-100:find(xDimMIP,1,'last')+100;
% imageRangeY = find(yDimMIP,1,'first')-100:find(yDimMIP,1,'last')+100;
% %Check to make sure imageRange is within image, otherwise adjust
% imageRangeX(imageRangeX < 1) = [];
% imageRangeX(imageRangeX > size(dicomStack,2)) = [];
% imageRangeY(imageRangeY < 1) = [];
% imageRangeY(imageRangeY > size(dicomStack,1)) = [];
% 
% figure
% for j=1:size(dicomStack,3)
%     ax1 = scrollsubplot(1,3,(j-1)*3+1); %%Need the scrollsubplot function in path
%     imagesc(dicomStack(imageRangeY,imageRangeX,j))
%     colormap(ax1, gray)
%     ax2 = scrollsubplot(1,3,(j-1)*3+2); %%Need the scrollsubplot function in path
%     imagesc(segVolume(imageRangeY,imageRangeX,j,1))
%     colormap(ax2,gray)
%     ax3 = scrollsubplot(1,3,(j-1)*3+3); %%Need the scrollsubplot function in path
%     imagesc(segVolumeDicoms(imageRangeY,imageRangeX,j,1))
%     colormap(ax3,gray)
% end
% 
% fprintf('QA Done. Press any key to continue\n')
% w = waitforbuttonpress
% fprint('Done\n')

