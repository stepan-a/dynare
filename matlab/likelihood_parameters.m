function [DynareModel, Q, H, fval, exit_flag, info, DLIK] = likelihood_parameters(xparam1,  BoundsInfo, DynareOptions, DynareModel, EstimatedParameters)

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

fval        = [];
exit_flag   = 1;
info        = zeros(4,1);
DLIK        = [];
Q           = [];
H           = [];

% Return, with endogenous penalty, if some parameters are smaller than the lower bound of the prior domain.
if isestimation(DynareOptions) && ~isequal(DynareOptions.mode_compute,1) && any(xparam1<BoundsInfo.lb)
    k = find(xparam1<BoundsInfo.lb);
    fval = Inf;
    exit_flag = 0;
    info(1) = 41;
    info(4)= sum((BoundsInfo.lb(k)-xparam1(k)).^2);
    if nargout>6 && DynareOptions.analytic_derivation
        DLIK=ones(length(xparam1),1);
    end
    return
end

% Return, with endogenous penalty, if some parameters are greater than the upper bound of the prior domain.
if isestimation(DynareOptions) && ~isequal(DynareOptions.mode_compute,1) && any(xparam1>BoundsInfo.ub)
    k = find(xparam1>BoundsInfo.ub);
    fval = Inf;
    exit_flag = 0;
    info(1) = 42;
    info(4)= sum((xparam1(k)-BoundsInfo.ub(k)).^2);
    if nargout>6 && DynareOptions.analytic_derivation
        DLIK=ones(length(xparam1),1);
    end
    return
end

% Get the covariance matrices for the structural innovations (Q) and the measurement error (H).
DynareModel = set_all_parameters(xparam1, EstimatedParameters, DynareModel);
Q = DynareModel.Sigma_e;
if ~DynareOptions.dsge_var
    H = DynareModel.H;
end

% Test if Q is positive definite.
if ~issquare(Q) || EstimatedParameters.ncx || isfield(EstimatedParameters,'calibrated_covariances')
    [Q_is_positive_definite, penalty] = ispd(Q(EstimatedParameters.Sigma_e_entries_to_check_for_positive_definiteness,EstimatedParameters.Sigma_e_entries_to_check_for_positive_definiteness));
    if ~Q_is_positive_definite
        fval = Inf;
        exit_flag = 0;
        info(1) = 43;
        info(4) = penalty;
        return
    end
    if isfield(EstimatedParameters,'calibrated_covariances')
        correct_flag=check_consistency_covariances(Q);
        if ~correct_flag
            penalty = sum(Q(EstimatedParameters.calibrated_covariances.position).^2);
            fval = Inf;
            exit_flag = 0;
            info(1) = 71;
            info(4) = penalty;
            return
        end
    end
end

% Test if H is positive definite.
if ~DynareOptions.dsge_var && (~issquare(H) || EstimatedParameters.ncn || isfield(EstimatedParameters,'calibrated_covariances_ME'))
    [H_is_positive_definite, penalty] = ispd(H(EstimatedParameters.H_entries_to_check_for_positive_definiteness,EstimatedParameters.H_entries_to_check_for_positive_definiteness));
    if ~H_is_positive_definite
        fval = Inf;
        info(1) = 44;
        info(4) = penalty;
        exit_flag = 0;
        return
    end
    if isfield(EstimatedParameters,'calibrated_covariances_ME')
        correct_flag=check_consistency_covariances(H);
        if ~correct_flag
            penalty = sum(H(EstimatedParameters.calibrated_covariances_ME.position).^2);
            fval = Inf;
            exit_flag = 0;
            info(1) = 72;
            info(4) = penalty;
            return
        end
    end
end