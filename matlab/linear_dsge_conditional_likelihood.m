function [fval, info, exit_flag, DLIK, Hess, SteadyState, trend_coeff, Model, DynareOptions, BayesInfo, DynareResults] = ...
        linear_dsge_conditional_likelihood(xparam1, DynareDataset, DatasetInfo, DynareOptions, DynareModel, EstimatedParameters, BayesInfo, BoundsInfo, DynareResults, derivatives_info)


% Initialization of the returned variables and others.
fval        = [];
SteadyState = [];
trend_coeff = [];
exit_flag   = 1;
info        = zeros(4,1);
DLIK        = [];
Hess       = [];


% Exit with error if analytical_derivation option is used.
if DynareOptions.analytic_derivation
    error('The analytic_derivation and conditional_likelihood are not compatible!')
end

%------------------------------------------------------------------------------
% 1. Get the structural parameters & define penalties
%------------------------------------------------------------------------------

[DynareModel, Q, H, fval, exit_flag, info] = likelihood_parameters(xparam1,  BoundsInfo, DynareOptions, DynareModel, EstimatedParameters);

if ~exit_flag
    return
end

Q_upper_chol = chol(Q);
iQ_upper_chol = chol(inv(Q));

% Return an error if the interface for measurement errors is used.
if ~isequal(H, zeros(size(H))) || EstimatedParameters.ncn || EstimatedParameters.ncx
    error('Option conditional_likelihood does not support declaration of measurement errors. You can specify the measurement errors in the model block directly by adding measurement equations.')
end

%------------------------------------------------------------------------------
% 2. call model setup & reduction program
%------------------------------------------------------------------------------

[T, R, SteadyState, exit_flag, info, DynareModel, DynareOptions, DynareResults] = reduced_form_model(DynareModel, DynareOptions, DynareResults);

if ~exit_flag
    return
end

% check endogenous prior restrictions
info = endogenous_prior_restrictions(T, R, DynareModel, DynareOptions, DynareResults);
if info(1)
    fval = Inf;
    info(4)=info(2);
    exit_flag = 0;
    return
end

% Define a vector of indices for the observed variables. Is this really usefull?...
BayesInfo.mf = BayesInfo.mf1;

% Define the constant vector of the measurement equation.
if DynareOptions.noconstant
    constant = zeros(DynareDataset.vobs, 1);
else
    if DynareOptions.loglinear
        constant = log(SteadyState(BayesInfo.mfys));
    else
        constant = SteadyState(BayesInfo.mfys);
    end
end

% Define the deterministic linear trend of the measurement equation.
if BayesInfo.with_trend
    [trend_addition, trend_coeff]=compute_trend_coefficients(DynareModel, DynareOptions, DynareDataset.vobs, DynareDataset.nobs);
    trend = repmat(constant, 1, DynareDataset.nobs)+trend_addition;
else
    trend_coeff = zeros(DynareDataset.vobs, 1);
    trend = repmat(constant, 1, DynareDataset.nobs);
end


% Return an error if some observations are missing.
if DatasetInfo.missing.state
    error('Option conditional_likelihood is not compatible with missing observations.')
end

% Get the selection matrix (vector of row indices for T and R)
Z = BayesInfo.mf;

% Get the number of observed variables.
pp = DynareDataset.vobs;

% Get the number of variables in the state equations (state variables plus observed variables).
mm = size(T, 1);

% Get the number of innovations.
rr = length(Q);

% Return an error if the number of shocks is not equal to the number of observations.
if ~isequal(pp, rr)
    error('With conditional_likelihood the number of innovations must be equal to the number of observed varilables!')
end

% Remove the trend.
Y = transpose(DynareDataset.data)-trend;

% Set state vector (deviation to steady state)
S = zeros(mm, 1);

%------------------------------------------------------------------------------
% 3. Evaluate the conditional likelihood
%------------------------------------------------------------------------------

Rtild = inv(R(Z,:));
const = -.5*rr*log(2*pi);
const = const + log(abs(det(Rtild))) + sum(log(diag(iQ_upper_chol)));

llik = zeros(size(Y, 2));

Ytild = Rtild*Y;
Ttild = Rtild*T(Z,:);

for t=1:size(Y, 2)
    epsilon = Ytild(:,t) - Ttild*S;
    upsilon = iQ_upper_chol*epsilon;
    S = T*S + R*epsilon;
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
    tmp = feval(DynareOptions.prior_restrictions.routine, DynareModel, DynareResults, DynareOptions, DynareDataset, DatasetInfo);
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