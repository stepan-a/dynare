function datatomfile (s, var_list, names)

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
%
% REMARKS 
% Only the first argument is mandatory. If only one input argument is provided, all the variables as defined in 
% M_.endo_names will be saved in the generated m file.

% Copyright (C) 2001-2018 Dynare Team
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

global M_ oo_

if nargin<2 || isempty(var_list)
    var_list = M_.endo_names(1:M_.orig_endo_nbr);
end

if nargin>2 && ~isempty(names)
    if ~isequal(length(var_list), length(names))
        error('datatomfile:: Second and third arguments must have the same number of rows (variables)!')
    end
else
    names = var_list;
    n = length(names);
end

data2mfile(s, var_list, names, M_, oo_);