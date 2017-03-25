function importfile(fileToRead)
%IMPORTFILE(FILETOREAD1)
%  Imports data from the specified file
%  FILETOREAD1:  file to read
 
%  Auto-generated by MATLAB on 25-Mar-2017 20:42:02
 
% Import the file
rawData1 = importdata(fileToRead);

% For some simple files (such as a CSV or JPEG files), IMPORTDATA might
% return a simple array.  If so, generate a structure so that the output
% matches that from the Import Wizard.
[~,name] = fileparts(fileToRead);
newData1.(matlab.lang.makeValidName(name)) = rawData1;

% Create new variables in the base workspace from those fields.
vars = fieldnames(newData1);
for i = 1:length(vars)
    assignin('base', vars{i}, newData1.(vars{i}));
end

