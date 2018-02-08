function [fval, info, exit_flag, DLIK, Hess, ys, trend_coeff, Model, DynareOptions, BayesInfo, DynareResults] = ...
        dsge_conditional_likelihood_2(xparam1, DynareDataset, DatasetInfo, DynareOptions, Model, EstimatedParameters, BayesInfo, BoundsInfo, DynareResults)

% Evaluates the conditional likelihood of a dsge model based on a second
% order local approximation, and return the posterior kernel.

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

% Initialization of the returned variables and others.
fval        = [];
SteadyState = [];
trend_coeff = [];
exit_flag   = 1;
info        = zeros(4,1);
DLIK        = [];
Hess        = [];

% Exit with error if analytical_derivation option is used.
if DynareOptions.analytic_derivation
    error('The analytic_derivation and conditional_likelihood are not compatible!')
end

% Issue an error if loglinear option is used.
if DynareOptions.loglinear
    error('It is not possible to use conditional_likelihood and order>1 with the option loglinear!')
end

% ----------------------------------------------------
%  1. Get the structural parameters & define penalties
% ----------------------------------------------------

% Return, with endogenous penalty, if some parameters are smaller than the lower bound of the prior domain.
if ~isequal(DynareOptions.mode_compute, 1) && any(xparam1<BoundsInfo.lb)
    k = find(xparam1<BoundsInfo.lb);
    fval = Inf;
    exit_flag = 0;
    info(1) = 41;
    info(4)= sum((BoundsInfo.lb(k)-xparam1(k)).^2);
    return
end

% Return, with endogenous penalty, if some parameters are greater than the upper bound of the prior domain.
if ~isequal(DynareOptions.mode_compute, 1) && any(xparam1>BoundsInfo.ub)
    k = find(xparam1>BoundsInfo.ub);
    fval = Inf;
    exit_flag = 0;
    info(1) = 42;
    info(4)= sum((xparam1(k)-BoundsInfo.ub(k)).^2);
    return
end

% Get the diagonal elements of the covariance matrices for the structural innovations (Q) and the measurement error (H).
Model = set_all_parameters(xparam1, EstimatedParameters, Model);

Q = Model.Sigma_e;
H = Model.H;

% Test if Q is positive definite.
if ~issquare(Q) || EstimatedParameters.ncx || isfield(EstimatedParameters,'calibrated_covariances')
    [Q_is_positive_definite, penalty] = ispd(Q(EstimatedParameters.Sigma_e_entries_to_check_for_positive_definiteness, EstimatedParameters.Sigma_e_entries_to_check_for_positive_definiteness));
    if ~Q_is_positive_definite
        fval = Inf;
        exit_flag = 0;
        info(1) = 43;
        info(4) = penalty;
        return
    end
    if isfield(EstimatedParameters, 'calibrated_covariances')
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

Q_upper_chol = chol(Q);
iQ_upper_chol = chol(inv(Q));

% Return an error if the interface for measurement errors is used.
if ~isequal(H, zeros(size(H))) || EstimatedParameters.ncn || EstimatedParameters.ncx
    error('conditional_likelihood does not support declaration of measurement errors. You can specify the measurement errors in the model block directly by adding measurement equations.')
end

% ----------------------------------------
%  2. call model setup & reduction program
% ----------------------------------------

% Linearize the model around the deterministic sdteadystate and extract the matrices of the state equation (T and R).
[T, R, SteadyState, info, Model, DynareOptions, DynareResults] = dynare_resolve(Model, DynareOptions, DynareResults,'restrict');

if info(1)
    if info(1) == 3 || info(1) == 4 || info(1) == 5 || info(1)==6 ||info(1) == 19 || ...
                info(1) == 20 || info(1) == 21 || info(1) == 23 || info(1) == 26 || ...
                info(1) == 81 || info(1) == 84 ||  info(1) == 85
        %meaningful second entry of output that can be used
        fval = Inf;
        info(4) = info(2);
        exit_flag = 0;
        return
    else
        fval = Inf;
        info(4) = 0.1;
        exit_flag = 0;
        return
    end
end

% Define the deterministic linear trend of the measurement equation.
if DynareOptions.noconstant
    constant = zeros(DynareDataset.vobs, 1);
else
    constant = SteadyState(BayesInfo.mfys);
end

% Define the deterministic linear trend of the measurement equation.
if BayesInfo.with_trend
    [trend_addition, trend_coeff] = compute_trend_coefficients(Model, DynareOptions, DynareDataset.vobs, DynareDataset.nobs);
    trend = repmat(constant, 1, DynareDataset.info.ntobs) + trend_addition;
else
    trend = repmat(constant, 1, DynareDataset.nobs);
end

% -------------------------------------------------
%  3. Get transition rules and transition equations
% -------------------------------------------------

dr = DynareResults.dr;
mf0 = BayesInfo.mf0;
mf1 = BayesInfo.mf1;
restrict_variables_idx  = dr.restrict_var_list;
observed_variables_idx  = restrict_variables_idx(mf1);
state_variables_idx     = restrict_variables_idx(mf0);
number_of_state_variables = length(mf0);
number_of_observed_variables = length(mf1);
number_of_structural_innovations = length(Q);

ghx  = dr.ghx(observed_variables_idx,:);
ghu  = dr.ghu(observed_variables_idx,:);
ghxx = dr.ghxx(observed_variables_idx,:);
ghuu = dr.ghuu(observed_variables_idx,:);
ghxu = dr.ghxu(observed_variables_idx,:);
steadystate = dr.ys(dr.order_var(observed_variables_idx));
constant = steadystate + .5*dr.ghs2(observed_variables_idx);

ghx_  = dr.ghx(state_variables_idx,:);
ghu_  = dr.ghu(state_variables_idx,:);
ghxx_ = dr.ghxx(state_variables_idx,:);
ghuu_ = dr.ghuu(state_variables_idx,:);
ghxu_ = dr.ghxu(state_variables_idx,:);
steadystate_ = dr.ys(dr.order_var(state_variables_idx));
constant_ = steadystate_ + .5*dr.ghs2(state_variables_idx);


% -------------------------
%  4. Likelihood evaluation
% -------------------------

% Get data and detrend.
Y = transpose(DynareDataset.data)-trend;

% Initialize the vector of logged densities.
llik = zeros(columns(Y));

% Initialize vector of state variables.
S = zeros(length(mf0), 1);

%  
for t = 1:columns(Y)
    epsilon = local_state_space_inversion_2(Y(:,t), S, ghx, ghu, constant, ghxx, ghuu, ghxu, DynareOptions.threads);
    upsilon = iQ_upper_chol*epsilon;
    S = local_state_space_iteration_2(S, epsilon(:,t), ghx_, ghu_, constant_, ghxx_, ghuu_, ghxu_, DynareOptions.threads);
    llik(t) = const - .5*(upsilon'*upsilon);
end

% Computes minus log-likelihood.
likelihood = -sum(llik(DynareOptions.presample+1:size(Y, 2)));


% ------------------------------------------------------------------------------
% 5. Adds prior if necessary
% ------------------------------------------------------------------------------

lnprior = priordens(xparam1, BayesInfo.pshape, BayesInfo.p6, BayesInfo.p7, BayesInfo.p3, BayesInfo.p4);

if DynareOptions.endogenous_prior==1
    [lnpriormom]  = endogenous_prior(Y, Pstar, BayesInfo, H);
    fval = (likelihood-lnprior-lnpriormom);
else
    fval = (likelihood-lnprior);
end

if DynareOptions.prior_restrictions.status
    tmp = feval(DynareOptions.prior_restrictions.routine, Model, DynareResults, DynareOptions, DynareDataset, DatasetInfo);
    fval = fval - tmp;
end

if isnan(fval)
    fval = Inf;
    info(1) = 47;
    info(4) = 0.1;
    exit_flag = 0;
    return
end

if imag(fval)~=0
    fval = Inf;
    info(1) = 48;
    info(4) = 0.1;
    exit_flag = 0;
    return
end