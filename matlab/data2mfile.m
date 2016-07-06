function data2mfile(s, var_list, names, M_, oo_)

% This optional command saves the simulation results in a text file. The name of each
% variable preceeds the corresponding results.
%
% INPUTS
%    s           [char]   Row char array, data file name.
%    var_list    [cell]   Cell of row char arrays, selected endogenous variables.
%    names       [cell]   Cell of row char arrays, alternative names for the endogenous variables in the data file.
%
% OUTPUTS
%    none
%
% SPECIAL REQUIREMENTS
%    none

% Copyright (C) 2018 Dynare Team
%
% This file is part of Dynare.
%
% Dynare is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
%
% Dynare is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with Dynare.  If not, see <http://www.gnu.org/licenses/>.

% Set default values for the second argument
if nargin<2 || isempty(var_list)
    var_list = M_.endo_names(1:M_.orig_endo_nbr);
end

% Set default values for the third argument or check consistency with the second argument
if nargin>2 && ~isempty(names)
    if ~isequal(length(var_list), length(names))
        error('Second and third arguments must have the same number of elements!')
    end
else
    names = var_list;
end

% Open the data file.
sm=[s,'.m'];
fid=fopen(sm,'w') ;

n = length(var_list);
ivar = zeros(n, 1);

% Get indices for the endogenous variables.
for i = 1:n
    i_tmp = strmatch(var_list{i}, M_.endo_names, 'exact');
    if isempty(i_tmp)
        error (['One of the specified variables does not exist']) ;
    else
        ivar(i) = i_tmp;
    end
end

% Save the selected data.
for i = 1:n
    fprintf(fid,[names{i}, ' = ['], '\n') ;
    fprintf(fid,'\n') ;
    fprintf(fid,'%15.8g\n', oo_.endo_simul(ivar(i),:)') ;
    fprintf(fid,'];\n') ;
    fprintf(fid,'\n') ;
end

% Close the data file.
fclose(fid);