function [T, R, SteadyState, exit_flag, info, DynareModel, DynareOptions, DynareResults, DLIK] = reduced_form_model(DynareModel, DynareOptions, DynareResults)

% Returns the reduced form model.

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


% Initialization of the returned arguments.
T = [];
R = [];
SteadyState = [];
exit_flag = 1;
info = [];

if nargout>8
    DLIK = [];
end


% Linearize the model around the deterministic steady state and extract the
% matrices of the state equation (T and R). Note that if order>1, the rest of
% the reduced form matrices are available in DynareResults.dr
[T, R, SteadyState, info, DynareModel, DynareOptions, DynareResults] = dynare_resolve(DynareModel, DynareOptions, DynareResults, 'restrict');

% Return, with endogenous penalty when possible, if dynare_resolve issues an error code (defined in resol).
if info(1)
    if info(1) == 3 || info(1) == 4 || info(1) == 5 || info(1)==6 ||info(1) == 19 ||...
                info(1) == 20 || info(1) == 21 || info(1) == 23 || info(1) == 26 || ...
                info(1) == 81 || info(1) == 84 ||  info(1) == 85 ||  info(1) == 86
        fval = Inf;
        info(4) = info(2);
        exit_flag = 0;
        if nargout>8 && DynareOptions.analytic_derivation
            DLIK=ones(length(xparam1),1);
        end
        return
    else
        fval = Inf;
        info(4) = 0.1;
        exit_flag = 0;
        if nargout>8 && DynareOptions.analytic_derivation
            DLIK=ones(length(xparam1),1);
        end
        return
    end
end