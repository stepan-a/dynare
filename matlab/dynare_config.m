function dynareroot = dynare_config(path_to_dynare, verbose)

% This function tests the existence of valid mex files (for qz
% decomposition, solution to sylvester equation and kronecker
% products...) and, if needed, add paths to the matlab versions
% of these routines.
% Also adds other directories to the path.
%
% INPUTS
%   none
%
% OUTPUTS
%   none
%
% SPECIAL REQUIREMENTS
%   none

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

if nargin && ~isempty(path_to_dynare)
    addpath(path_to_dynare);
end

dynareroot = strrep(which('dynare'),'dynare.m','');

origin = pwd();
cd([dynareroot '/..'])

if ~nargin || nargin==1
    verbose = 1;
end

p = {'/distributions/' ; ...
     '/kalman/' ; ...
     '/kalman/likelihood' ; ...
     '/AIM/' ; ...
     '/partial_information/' ; ...
     '/perfect-foresight-models/' ; ...
     '/ms-sbvar/' ; ...
     '/ms-sbvar/identification/' ; ...
     '/../contrib/ms-sbvar/TZcode/MatlabFiles/' ; ...
     '/../contrib/jsonlab/' ; ...
     '/parallel/' ; ...
     '/particles/src' ; ...
     '/gsa/' ; ...
     '/ep/' ; ...
     '/backward/' ; ...
     '/convergence_diagnostics/' ; ... 
     '/cli/' ; ...
     '/lmmcp/' ; ...
     '/optimization/' ; ...
     '/modules/dates/src/' ; ...
     '/modules/dseries/src/' ; ...
     '/utilities/doc/' ; ...
     '/utilities/tests/src/' ; ...
     '/utilities/dataset/' ; ...
     '/utilities/general/' ; ...
     '/utilities/graphics/' ; ...
     '/modules/reporting/src/' ; ...
     '/modules/reporting/macros/'; ...
     '/reduced_form_models/local_state_space_inversion_2/'};

% For functions that exist only under some Octave versions
% or some MATLAB versions, and for which we provide some replacement functions

if ~isoctave
    % Replacements for rows(), columns() and issquare() (inexistent under MATLAB)
    p{end+1} = '/missing/rows_columns';
    p{end+1} = '/missing/issquare';
    % Replacement for vec() (inexistent under MATLAB)
    p{end+1} = '/missing/vec';
    if ~user_has_matlab_license('statistics_toolbox')
        % Replacements for functions of the stats toolbox
        p{end+1} = '/missing/stats/';
    end
end

% ordeig() doesn't exist in Octave
if isoctave
    p{end+1} = '/missing/ordeig';
end

if isoctave && ~isequal(supported_octave_version(), version())
    skipline()
    warning(['This version of Octave is not supported. Consider installing ' ...
             'version %s of Octave,\notherwise m files will be used instead ' ...
             'of precompiled mex files and some features, like solution\n' ...
             'of models approximated at third order, will not be available.'], supported_octave_version())
    skipline()
end

% ilu is missing in Octave < 4.0
if isoctave && octave_ver_less_than('4.0')
    p{end+1} = '/missing/ilu';
end

% corrcoef with two outputs is missing in Octave (ticket #796)
if isoctave && ~user_has_octave_forge_package('nan')
    p{end+1} = '/missing/corrcoef';
end

% nanmean is in Octave Forge Statistics package and in MATLAB Statistics
% toolbox
if (isoctave && ~user_has_octave_forge_package('statistics')) ...
        || (~isoctave && ~user_has_matlab_license('statistics_toolbox'))
    p{end+1} = '/missing/nanmean';
end

% Check if struct2array is available.
if ~exist('struct2array')
    p{end+1} = '/missing/struct2array';
end

% isfile is missing in Octave and Matlab<R2017b
if isoctave || matlab_ver_less_than('9.3')
    p{end+1} = '/missing/isfile';
end

P = cellfun(@(c)[dynareroot(1:end-1) c], p, 'uni',false);

% Get mex files folder(s)
mexpaths = add_path_to_mex_files(dynareroot, false);

% Add mex files folder(s)
P(end+1:end+length(mexpaths)) = mexpaths;

% Set matlab's path
addpath(P{:});

% Set mex routine names
mex_status = cell(1,3);
mex_status(1,1) = {'mjdgges'};
mex_status(1,2) = {'qz'};
mex_status(1,3) = {'Generalized QZ'};
mex_status(2,1) = {'gensylv'};
mex_status(2,2) = {'gensylv'};
mex_status(2,3) = {'Sylvester equation solution'};
mex_status(3,1) = {'A_times_B_kronecker_C'};
mex_status(3,2) = {'kronecker'};
mex_status(3,3) = {'Kronecker products'};
mex_status(4,1) = {'sparse_hessian_times_B_kronecker_C'};
mex_status(4,2) = {'kronecker'};
mex_status(4,3) = {'Sparse kronecker products'};
mex_status(5,1) = {'local_state_space_iteration_2'};
mex_status(5,2) = {'reduced_form_models/local_state_space_iteration_2'};
mex_status(5,3) = {'Local state space iteration (second order)'};
number_of_mex_files = size(mex_status,1);

% Remove some directories from matlab's path. This is necessary if the user has
% added dynare_v4/matlab with the subfolders. Matlab has to ignore these
% subfolders if valid mex files exist.
matlab_path = path;
for i=1:number_of_mex_files
    test = strfind(matlab_path,[dynareroot mex_status{i,2}]);
    action = length(test);
    if action
        rmpath([dynareroot mex_status{i,2}]);
        matlab_path = path;
    end
end

% Test if valid mex files are available, if a mex file is not available
% a matlab version of the routine is included in the path.
if verbose
    skipline()
    disp('Configuring Dynare ...')
end

for i=1:number_of_mex_files
    test = (exist(mex_status{i,1},'file') == 3);
    if ~test
        addpath([dynareroot mex_status{i,2}]);
        message = '[m]   ';
    else
        message = '[mex] ';
    end
    if verbose
        disp([ message mex_status{i,3} '.' ])
    end
end

% Test if bytecode DLL is present
if exist('bytecode', 'file') == 3
    message = '[mex] ';
else
    message = '[no]  ';
end
if verbose
    disp([ message 'Bytecode evaluation.' ])
end

% Test if k-order perturbation DLL is present
if exist('k_order_perturbation', 'file') == 3
    message = '[mex] ';
else
    message = '[no]  ';
end
if verbose
    disp([ message 'k-order perturbation solver.' ])
end

% Test if dynare_simul_ DLL is present
if exist('dynare_simul_', 'file') == 3
    message = '[mex] ';
else
    message = '[no]  ';
end
if verbose
    disp([ message 'k-order solution simulation.' ])
end

% Test if qmc_sequence DLL is present
if exist('qmc_sequence', 'file') == 3
    message = '[mex] ';
else
    message = '[no]  ';
end
if verbose
    disp([ message 'Quasi Monte-Carlo sequence (Sobol).' ])
end

% Test if MS-SBVAR DLL is present
if exist('ms_sbvar_command_line', 'file') == 3
    message = '[mex] ';
else
    message = '[no]  ';
end
if verbose
    disp([ message 'Markov Switching SBVAR.' ])
    skipline()
end

% Initialization of the dates and dseries classes (recursive).
initialize_dseries_toolbox;

cd(origin);