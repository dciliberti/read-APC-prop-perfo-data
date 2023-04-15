function position = headFun(cellArray,stringa)
% Expects a 1×17 cell array and a character array as input
% 
% Input example:
%   Columns 1 through 8
%     {'RPM'}    {'V_kph'}    {'V_mph'}    {'J'}    {'eff'}    {'CT'}    {'CP'}    {'P_hp'}
% 
%   Columns 9 through 15
%     {'Q_lbfin'}    {'T_lbf'}    {'P_W'}    {'Q_Nm'}    {'T_N'}    {'T_P_g_over_W'}    {'Mach'}
% 
%   Columns 16 through 17
%     {'Reynolds'}    {'FOM'}

% Compare input character vector with all character vector in cell array
% (not case sensitive) and locate its position
confronto = zeros(size(cellArray,2),1);
for idx = 1:size(cellArray,2)
    confronto(idx) = strcmpi(stringa,cellArray{idx});
end

if any(confronto)
    position = find(confronto);
else
    error(['ERROR: The character array ', stringa, ' was not found as in cell array.'])
%     disp('Allowable names:')
%     disp(char(cellArray))
end
