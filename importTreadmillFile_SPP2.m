function treadmillData = importTreadmillFile(filename, dataLines)
%IMPORTFILE Import data from a text file
%  CONDITION1TREADMILL0010001 = IMPORTFILE(FILENAME) reads data from
%  text file FILENAME for the default selection.  Returns the data as a
%  table.
%
%  CONDITION1TREADMILL0010001 = IMPORTFILE(FILE, DATALINES) reads data
%  for the specified row interval(s) of text file FILENAME. Specify
%  DATALINES as a positive scalar integer or a N-by-2 array of positive
%  scalar integers for dis-contiguous row intervals.
%
%  Example:
%  condition1treadmill0010001 = importfile("/Users/hjhuang/Dropbox (Personal)/_MANUSCRIPTS/_Cesar/SP controller testing data/condition1_treadmill0010001.txt", [2, Inf]);
%
%  See also READTABLE.
%
% Auto-generated by MATLAB on 19-Mar-2021 00:20:21

%% Input handling

% If dataLines is not specified, define defaults
if nargin < 2
    dataLines = [2, Inf];
end

%% Setup the Import Options and import the data
opts = delimitedTextImportOptions("NumVariables", 18);

% Specify range and delimiter
opts.DataLines = dataLines;
opts.Delimiter = "\t";

% Specify column names and types
opts.VariableNames = ["Time", "t", "stride", "SpeedTheory1", "SpeedTheory2", "x1", "x", "SwayTheory", "COMz", "LFTa", "RFTa", "LKNEa", "RKNEa", "LHIPa", "RHIPa", "SpeedActual1", "SpeedActual2", "SwayActual"];
opts.VariableTypes = ["double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double"];

% Specify file level properties
opts.ExtraColumnsRule = "ignore";
opts.EmptyLineRule = "read";

% Import the data
treadmillData = readtable(filename, opts);

end