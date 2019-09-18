function [stringnames] = detectfiles(folderdir,detectstring)

stringnames = {};
counter = 1;
folderinfo = dir(folderdir);
for i = 1:length(folderinfo)
    if contains(folderinfo(i).name,detectstring) && folderinfo(i).name(1)~='.'
        stringnames{counter} = folderinfo(i).name;
        counter = counter+1;
    end
end

end