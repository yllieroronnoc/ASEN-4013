% savefig.m
% Saves current figure as .fig and .tiff files, in the current directory,
% using the chart figure.title as the file name
% Travis Sippel
% tsippel@purdue.edu
% Oct 30, 2009

% if tif='y' or 'Y', a tif image will be created for each of the figure
% numbers in 'handles'
% if fig='y' or 'Y', a matlab .fig file will be created for each of the
% figure numbers in the array 'handles'

% Example: if you want to save figure numbers 5, 7, and 14 to tif file only
% in the current matlab directory, use the following commands:
% savefig(5,'file1','y','n');
% savefig(7,'file2','y','n');
% savefig(14,'file3','y','n');

function savefig(handles, filename, tif, fig)

for i = 1:length(handles);
    if tif == 1 || tif == 'y' || tif == 'Y';
        saveas(handles(i),[filename, '.tif']);
    end;
    if fig == 1 || fig == 'y' || fig == 'Y'
        saveas(handles(i),[filename, '.fig']);
    end;
end;

