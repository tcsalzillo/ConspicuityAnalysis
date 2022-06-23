function main
%Run this function to generate conspicuity measurements in batch. Ensure 
%that importImage, importSeg, measureConspicuity, optimizeCentroid, and 
%combineConspicuityData are also located in this folder. 
%
%Folder heirarchy should be:
%   This folder (includes .m files and "Data" folder)
%       Data folder with patient MRN subfolders
%           MRN folder with image subfolders
%               Image folder with "dicoms" subfolder and one subfolder
%               for each structure
%                   The "dicoms" folder should contain the image dicoms
%                   The structure folders should contain the structures and
%                   corresponding expanded/contracted structures, all in
%                   .nii.gz format
%
%Ensure that the following mrnList, MRList, structList, and obList are
%correct for your data.

flag = 0;
imageImport = 0; %Set to 1 to import image
segmentationImport = 0; %Set to 1 to import segmentations
conspicuityTool = 0; %Set to 1 to measure conspicuity
writeConspicuity = 0; %Set to 1 to write conspicuity measurements to CSV files
gatherConspictuityData = 1; %Set to 1 to combine individual conspicuity CSV files 
                            %into 1 master CSV file. These .csv files need
                            %to be copied from the individual structure
                            %folders to the ConspicuityMeasurements folder
                            %according to the following heirarchy:
                            %"ConspicuityMeasurements" -> MRN Folder -> 
                            %Image Folder -> all structure .csv files. 
                            %Thus, this should only be set to 1 after the 
                            %conspicuity .csv files have been generated 
                            %(hint: turn the other values to 0 when you 
                            %turn this to 1 and vice versa)

%Set master directory
superDir = sprintf('%s%s',pwd,'\Data\'); %calls current directory- change if needed

mrnList = {'Pat1'}; %separate with semicolons
numMRN = size(mrnList,1);

MRList = {'Image1'};
numMR = size(MRList,1);

structList = {'Structure1'};
numStruct = size(structList,1);

obList = {'Ob1'};
numObservers = size(obList,1);

conspicuityDataDir = sprintf('%s%s',superDir,'ConspicuityMeasurements\');
compiledData = {'Label','Median Conspicuity','Mean Conspicuity', 'Mean Conpsicuity (90th%)'};
allConspData = {'Label','Image','All Conspicuity Values'};
missingMeasurements = {};

%PixShrink refers to number of pixels to shrink the structure seg to produce
%the innerSegShell. By default, this is set to 2, but for small structures
%such as lymph nodes, you may only want to shrink by 1 pixel so that there
%are sufficient pixels in the innerSegShell. PixGrow can also be adjusted
%to 1 or 2, depending on how far away from the segmentation you wish to
%measure conspicuity. By default, this is set to 2 to allow for ~1 mm
%tolerance of segmentations
pixShrink = 2;
pixGrow = 2;

for w = 1:numMRN %Loop through MRN
    for x = 1:numMR %Loop through images
        mrn = char(mrnList(w,1));
        MRName = char(MRList(x,1));
        
        %Check for missing MRN and MR Folders
        if imageImport == 1
            if ~any(size(dir(sprintf('%s%s',superDir,mrn)),1))
                fprintf('%s%s%s\n','MRN: ',mrn,' folder not present') 
                continue
            end
            if ~any(size((dir(sprintf('%s%s%s%s',superDir,mrn,'\',MRName))),1))
                fprintf('%s%s\n',MRName,' folder not present') 
                continue
            end
            
            dirName = sprintf('%s%s%s%s', superDir,mrn,'\',MRName);
            
            %Import Image
            [dicomStack,flag] = importImage(dirName); %make sure this .m file is in path
            if flag == 1
                fprintf('%s%s\n','No image data present in MR:',MRName)
                continue
            end
        end
        
        for y = 1:numStruct %Loop through structures
            for z = 1:numObservers %Loop through observers                 
                structure = char(structList(y,1));
                observerName = char(obList(z,1)); 
                fprintf('\n%s%s%s%s%s%s%s%s\n','MRN: ',mrn,'  MR: ',MRName,'  Structure: ',structure,'  Observer: ',observerName);
                
                %Import Structure Segmentations
                if ((imageImport == 1) && (segmentationImport == 1))
                    [segVolume,segShell,segVolumeDicoms,segShellDicoms,flag] = importSeg(dirName,structure,observerName,dicomStack,pixShrink,pixGrow);
%                     if isempty(segVolume)
%                         fprintf('%s%s%s%s\n','No structure segmentation data present in MR:',MRNum,'  Observer:',caseNum)
%                         continue
%                     end 
                    if flag == 1
                        fprintf('Missing sub-structure(s) (check filenames)\n')
                        continue
                    end                   
                end
     
                %Run Conspicuity Measurement
                if ((imageImport == 1) && (segmentationImport == 1) && (conspicuityTool == 1))
                    [conspicuityData,flag] = measureConspicuity(dicomStack, segVolume,segShell,segVolumeDicoms,segShellDicoms);

                    tableName = sprintf('%s%s%s%s%s%s%s%s',dirName,'\',structure,'\',MRName,structure,observerName,'ConspicuityData.csv');
                    if flag == 1
                        fprintf('Segmentation empty in %s\n',tableName)
                        continue
                    elseif flag == 2
                        fprintf('Conspicuity Matrix empty in %s\n',tableName)
                        continue
                    elseif writeConspicuity == 1
                        fprintf('Mean of inner 90th percentile conspicuity: %d \n',str2double(conspicuityData{2,4}))
                        if exist(tableName,'file')==2
                            delete(tableName);
                        end
                        conspicuityTable = table(conspicuityData);
                        writetable(conspicuityTable,tableName,'Delimiter',',','QuoteStrings',true);
                        fprintf('Conspicuity Data written to %s\n',tableName)
                    else
                        fprintf('Mean of inner 90th percentile conspicuity: %d \n',str2double(conspicuityData{2,4}))
                    end
                end    
                
                %Run Conspicuity Data Compilation
                if gatherConspictuityData == 1
                    conspicuityDirName = sprintf('%s%s%s%s', conspicuityDataDir,mrn,'\',MRName);
                    tableName = sprintf('%s%s%s%s',MRName,structure,observerName,'ConspicuityData.csv');  
                    
                    [dataEntry,allConspEntry,flag] = combineConpsicuityData(conspicuityDirName,mrn,MRName,tableName);
                    if flag == 1
                        fprintf('%s file not present\n',tableName)
                        missingMeasurements = [missingMeasurements;tableName];
%                         pause(0.3)
                        continue
                    else
                        compiledData = [compiledData;dataEntry]; 
                        allConspData = [allConspData;allConspEntry];
                    end
                end            
            end      
        end
    end
end

if gatherConspictuityData == 1
    compiledTable = table(compiledData); %Non-organized conspicuity data in individual rows
    compiledTableName = sprintf('%s%s',conspicuityDataDir,'CompiledConspicuityData.csv');
    if exist(compiledTableName,'file')==2
        delete(compiledTableName);
    end
    writetable(compiledTable,compiledTableName,'Delimiter',',','QuoteStrings',true);
    
    allConspTable = table(allConspData); %Non-organized conspicuity data in individual rows
    allConspTableName = sprintf('%s%s',conspicuityDataDir,'allConspicuityData.csv');
    if exist(compiledTableName,'file')==2
        delete(compiledTableName);
    end
    writetable(allConspTable,allConspTableName,'Delimiter',',','QuoteStrings',true);
    
    %Generate CSV to display missing conspicuity measurements
    missingTable = table(missingMeasurements);
    missingTableName = sprintf('%s%s',conspicuityDataDir,'MissingConspicuityData.csv');
    if exist(missingTableName,'file')==2
        delete(missingTableName);
    end
    writetable(missingTable,missingTableName,'Delimiter',',','QuoteStrings',true);
end