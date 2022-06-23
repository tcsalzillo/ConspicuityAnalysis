function[dicomStack,flag] = importImage(dirName)
%This function imports DICOM images located in "\dicoms" subdirectory of
%dirName and outputs a matrix (dicomStack) of the grayscale values

dicomStack = [];
flag = 0;
fprintf('Importing Image from %s\n',dirName);

dicoms = dir(sprintf('%s%s',dirName,'\dicoms'));
numIm = size(dicoms,1)-2; %For some reason, first 2 files, which appear hidden,
                            % aren't images;
if ~any(size(dir([sprintf('%s%s',dirName,'\dicoms') '/*.dcm' ]),1))
%     fprintf('No DICOMs present')
    flag = 1;
    return
end

dicomName = sprintf('%s%s%s',dicoms(end).folder,'\',dicoms(end).name);
temp1 = dicomread(dicomName);
dicomStack = zeros(size(temp1,1),size(temp1,2),numIm);
%Creates dicomStack size using last DICOM image as template

for j = 1:numIm
    dicomName = sprintf('%s%s%s',dicoms(j+2).folder,'\',dicoms(j+2).name);
    %Again, first two files aren't images
    tempDicomStack(:,:,j) = dicomread(dicomName);
end

%Below might change between scanners or image types- verify by uncommenting
%the scrollsubplot figure and checking the image stack manually
dicomStack = circshift(tempDicomStack,floor(numIm/2),3); %MR images started halfway through the stack