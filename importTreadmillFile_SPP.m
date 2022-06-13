function s20nopertrubtreadmill0001 = importTreadmillFile_SPP(filename, dataLines)
%IMPORTFILE Import data from a text file
%  S20NOPERTRUBTREADMILL0001 = IMPORTFILE(FILENAME) reads data from text
%  file FILENAME for the default selection.  Returns the data as a table.
%
%  S20NOPERTRUBTREADMILL0001 = IMPORTFILE(FILE, DATALINES) reads data
%  for the specified row interval(s) of text file FILENAME. Specify
%  DATALINES as a positive scalar integer or a N-by-2 array of positive
%  scalar integers for dis-contiguous row intervals.
%
%  Example:
%  s20nopertrubtreadmill0001 = importfile("F:\SPP\Helen\s20_no_pertrub_treadmill0001.txt", [2, Inf]);
%
%  See also READTABLE.
%
% Auto-generated by MATLAB on 08-Oct-2021 17:46:32

%% Input handling

% If dataLines is not specified, define defaults
if nargin < 2
    dataLines = [2, Inf];
end

%% Set up the Import Options and import the data
opts = delimitedTextImportOptions("NumVariables", 20);

% Specify range and delimiter
opts.DataLines = dataLines;
opts.Delimiter = "\t";

% Specify column names and types
opts.VariableNames = ["Time", "t", "stride", "SpeedTheory1", "SpeedTheory2", "x1", "x", "SwayTheory", "COMz", "LFTa", "RFTa", "LKNEa", "RKNEa", "LHIPa", "RHIPa", "SpeedActual1", "SpeedActual2", "SwayActual", "LeftBeltDistance", "RightBeltDistance"];
opts.VariableTypes = ["double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double"];

% Specify file level properties
opts.ExtraColumnsRule = "ignore";
opts.EmptyLineRule = "read";

% Import the data
s20nopertrubtreadmill0001 = readtable(filename, opts);

end