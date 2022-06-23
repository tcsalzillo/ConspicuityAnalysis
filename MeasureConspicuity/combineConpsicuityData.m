function[dataEntry,allConspEntry,flag] = combineConpsicuityData(conspicuityDirName,mrn,MRName,tableName)
flag = 0;
dataEntry = {};
allConspEntry = {};

if ~any(size(dir(sprintf('%s%s%s',conspicuityDirName,'\',tableName)),1))
    flag = 1;
    return
end

imageNum = tableName(1:size(MRName,2));
               
oldDir = cd(conspicuityDirName);
conspicuityData = readtable(tableName);

minConspicuity = conspicuityData{1,2};
maxConspicuity = conspicuityData{1,4};
meanConspicuity = conspicuityData{2,2};
meanConspicuity90 = conspicuityData{2,4};
rawConspicuity = conspicuityData{4:end,4};
medianConspicuity = median(rawConspicuity);

dataLabel = sprintf('%s%s',mrn,tableName);
dataEntry = {dataLabel,num2str(medianConspicuity,'%.4f')...
    num2str(meanConspicuity,'%.4f'),num2str(meanConspicuity90,'%.4f')};
allConspEntry = {dataLabel,imageNum,num2str(rawConspicuity')};

cd(oldDir);