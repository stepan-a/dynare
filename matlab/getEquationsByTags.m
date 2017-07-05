function [lhs, rhs, linenum] = getEquationsByTags(jsonmodel, tagname, tagvalue)
%function [lhs, rhs] = getEquationByTag(jsonmodel, tag)
% Return the lhs, rhs of an equation and the line it was defined
% on given its tag
%
% INPUTS
%   jsonmodel
%   tagname
%   tagvalue
%
% OUTPUTS
%   lhs
%   rhs
%   linenum
%
% SPECIAL REQUIREMENTS
%   none

% Copyright (C) 2017 Dynare Team
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

if ~iscell(tagname)
    tagname = {tagname};
end
if ~iscell(tagvalue)
    tagvalue = {tagvalue};
end

assert(length(tagname) == length(tagvalue));

if length(tagname) == 1
    lhs = '';
    rhs = '';
    linenum = -1;
else
    lhs = cell(1, length(tagname));
    rhs = cell(1, length(tagname));
    linenum = cell(1, length(tagname));
end

for i=1:length(jsonmodel)
    for j = 1:length(tagname)
        if isfield(jsonmodel{i}, 'tags') && ...
                isfield(jsonmodel{i}.tags, tagname{j}) && ...
                strcmp(jsonmodel{i}.tags.(tagname{j}), tagvalue{j})
            if length(tagname) == 1
                lhs = jsonmodel{i}.lhs;
                rhs = jsonmodel{i}.rhs;
                linenum = jsonmodel{i}.line;
                return
            else
                lhs{j} = jsonmodel{i}.lhs;
                rhs{j} = jsonmodel{i}.rhs;
                linenum{j} = jsonmodel{i}.line;
                if ~any(cellfun(@isempty, lhs))
                    return
                end
                break
            end
        end
    end
end
end