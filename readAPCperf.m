function apcPropPerfData = readAPCperf(filename, startRow, endRow)
%READAPCPERF Import APC propeller performance data from file.
% APC performance data files usually report data as tables at different
% RPM values. All data are collected in a table, so that the variables
% can be extracted when needed with dot notation. Therefore, this 
% function distinguishes each RPM value, separating the different
% datasets and converting the table into a cell array.
% Then, each dataset is moved into a cell of the cell array. Finally, 
% it converts each cell of cell array into a table and keeps all the 
% tables into the cell array adding columns for RPM and velocity in km/h.
% That is, the output of the function is a cell array containing several
% tables, where each table contains propeller performance data at a
% specific RPM value.
%
%   APCPROPPERFDATA = READAPCPERF(FILENAME)
%   Reads data from text file READAPCPERF for the default selection.
%
%   APCPROPPERFDATA = READAPCPERF(FILENAME, STARTROW, ENDROW)
%   Reads data from rows STARTROW through ENDROW of text file FILENAME.
%
% Example:
%   APCPROPPERFDATA = READAPCPERF('PER3_19x16.dat', 1, 463);
%
%    See also TEXTSCAN.

% Auto-generated by MATLAB on 2023/03/13 12:30:13
% Edited by Danilo Ciliberti on 2023/03/14

%% Initialize variables.
if nargin<=2
    startRow = 1;
    endRow = inf;
end

%% Read columns of data as text:
% For more information, see the TEXTSCAN documentation.
formatSpec = '%12s%12s%12s%12s%12s%12s%12s%12s%12s%12s%12s%12s%12s%12s%10s%s%[^\n\r]';

%% Open the text file.
fileID = fopen(filename,'r');

%% Read columns of data according to the format.
% This call is based on the structure of the file used to generate this code. If an error occurs for a different file, try regenerating the code from the Import Tool.
dataArray = textscan(fileID, formatSpec, endRow(1)-startRow(1)+1, 'Delimiter', '', 'WhiteSpace', '', 'TextType', 'string', 'HeaderLines', startRow(1)-1, 'ReturnOnError', false, 'EndOfLine', '\r\n');
for block=2:length(startRow)
    frewind(fileID);
    dataArrayBlock = textscan(fileID, formatSpec, endRow(block)-startRow(block)+1, 'Delimiter', '', 'WhiteSpace', '', 'TextType', 'string', 'HeaderLines', startRow(block)-1, 'ReturnOnError', false, 'EndOfLine', '\r\n');
    for col=1:length(dataArray)
        dataArray{col} = [dataArray{col};dataArrayBlock{col}];
    end
end

%% Close the text file.
fclose(fileID);

%% Convert the contents of columns containing numeric text to numbers.
% Replace non-numeric text with NaN.
raw = repmat({''},length(dataArray{1}),length(dataArray)-1);
for col=1:length(dataArray)-1
    raw(1:length(dataArray{col}),col) = mat2cell(dataArray{col}, ones(length(dataArray{col}), 1));
end
numericData = NaN(size(dataArray{1},1),size(dataArray,2));

for col=[1,2,3,4,5,6,7,8,9,10,11,12,13,14,15]
    % Converts text in the input cell array to numbers. Replaced non-numeric text with NaN.
    rawData = dataArray{col};
    for row=1:size(rawData, 1)
        % Create a regular expression to detect and remove non-numeric prefixes and suffixes.
        regexstr = '(?<prefix>.*?)(?<numbers>([-]*(\d+[\,]*)+[\.]{0,1}\d*[eEdD]{0,1}[-+]*\d*[i]{0,1})|([-]*(\d+[\,]*)*[\.]{1,1}\d+[eEdD]{0,1}[-+]*\d*[i]{0,1}))(?<suffix>.*)';
        try
            result = regexp(rawData(row), regexstr, 'names');
            numbers = result.numbers;

            % Detected commas in non-thousand locations.
            invalidThousandsSeparator = false;
            if numbers.contains(',')
                thousandsRegExp = '^[-/+]*\d+?(\,\d{3})*\.{0,1}\d*$';
                if isempty(regexp(numbers, thousandsRegExp, 'once'))
                    numbers = NaN;
                    invalidThousandsSeparator = true;
                end
            end
            % Convert numeric text to numbers.
            if ~invalidThousandsSeparator
                numbers = textscan(char(strrep(numbers, ',', '')), '%f');
                numericData(row, col) = numbers{1};
                raw{row, col} = numbers{1};
            end
        catch
            raw{row, col} = rawData{row};
        end
    end
end


%% Split data into numeric and string columns.
rawNumericColumns = raw(:, [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15]);
rawStringColumns = string(raw(:, 16));


%% Replace non-numeric cells with NaN
R = cellfun(@(x) ~isnumeric(x) && ~islogical(x),rawNumericColumns); % Find non-numeric cells
rawNumericColumns(R) = {NaN}; % Replace non-numeric cells

%% Make sure any text containing <undefined> is properly converted to an <undefined> categorical
idx = (rawStringColumns(:, 1) == "<undefined>");
rawStringColumns(idx, 1) = "";

%% Create output variable
apcPropPerfData = table;
apcPropPerfData.V_mph = cell2mat(rawNumericColumns(:, 1));
apcPropPerfData.J = cell2mat(rawNumericColumns(:, 2));
apcPropPerfData.eff = cell2mat(rawNumericColumns(:, 3));
apcPropPerfData.CT = cell2mat(rawNumericColumns(:, 4));
apcPropPerfData.CP = cell2mat(rawNumericColumns(:, 5));
apcPropPerfData.P_hp = cell2mat(rawNumericColumns(:, 6));
apcPropPerfData.Q_lbfin = cell2mat(rawNumericColumns(:, 7));
apcPropPerfData.T_lbf = cell2mat(rawNumericColumns(:, 8));
apcPropPerfData.P_W = cell2mat(rawNumericColumns(:, 9));
apcPropPerfData.Q_Nm = cell2mat(rawNumericColumns(:, 10));
apcPropPerfData.T_N = cell2mat(rawNumericColumns(:, 11));
apcPropPerfData.T_P_g_over_W = cell2mat(rawNumericColumns(:, 12));
apcPropPerfData.Mach = cell2mat(rawNumericColumns(:, 13));
apcPropPerfData.Reynolds = cell2mat(rawNumericColumns(:, 14));
apcPropPerfData.FOM = cell2mat(rawNumericColumns(:, 15));
apcPropPerfData.VarName16 = categorical(rawStringColumns(:, 1));

%% Edit table
apcPropPerfData.VarName16 = []; % remove last empty column
apcPropPerfData(1:16,:) = []; % remove first 16 rows
[nRows, ~] = size(apcPropPerfData);

% Read table row by row and separate data with RPM
dataset = 0; % counter
perfo = cell(1,99); % empty cell array
for row = 1:nRows
    if any(isnan(apcPropPerfData{row,:})) % determine if row contains NaN
        if ~isnan(apcPropPerfData{row,3}) % if the element is numeric that is the RPM value for the dataset
            dataset = dataset + 1; % increase dataset counter
            rpm(dataset) = apcPropPerfData{row,3}; % keep in memory RPM for the dataset
        end
    else % dataset value should be at least 1
        perfo{:,dataset} = [perfo{:,dataset}; num2cell(apcPropPerfData{row,:})]; % populate cell array with numeric data
    end
end
perfo = perfo(~cellfun('isempty',perfo)); % remove empty cell arrays

% Convert each cell of cell array into a table then keep all tables into a 
% cell array and add columns for velocity in km/h and RPM
header = apcPropPerfData.Properties.VariableNames;
clear apcPropPerfData % we need to clear original variable name
for idx = 1:dataset
    apcPropPerfData{idx} = cell2table(perfo{idx});
    apcPropPerfData{idx}.Properties.VariableNames = header;
    apcPropPerfData{idx}.V_kph = apcPropPerfData{idx}.V_mph * 1.60934;
    apcPropPerfData{idx}.RPM = ones(size(apcPropPerfData{idx},1),1) * rpm(idx);
    apcPropPerfData{idx} = [apcPropPerfData{idx}(:,end), apcPropPerfData{idx}(:,end-1), apcPropPerfData{idx}(:,1:end-2)];
end

end