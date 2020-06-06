%Reading in files from a directory
Files=dir('*.*');
for k=1:length(Files)
   FileNames=Files(k).name
end