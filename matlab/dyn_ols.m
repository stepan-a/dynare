function dyn_ols(ds, varargin)
% function dyn_ols(ds, varargin)
% Run OLS on chosen model equations; unlike olseqs, allow for time t
% endogenous variables on LHS
%
% INPUTS
%   ds      [dseries]    data
%
% OUTPUTS
%   none
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

global M_ oo_

jsonfile = [M_.fname '_original.json'];
if exist(jsonfile, 'file') ~= 2
    error('Could not find %s! Please use the json option (See the Dynare invocation section in the reference manual).', jsonfile);
end

%% Get Equation(s)
jsonmodel = loadjson(jsonfile);
jsonmodel = jsonmodel.model;
[lhs, rhs, lineno] = getEquationsByTags(jsonmodel, 'name', varargin{:});

%% Estimation
regexpr1 = ...
    ['(diff\(\w+(\(\W?\w+\))?\))\*$' ...
    '|' '\((\w+(\(\W?\w+\))?(\W?\w+(\(\W?\w+\))?)*)\)\*$' ...
    '|' '(\w+(\(\W?\w+\))?)\*$' ...
    ];

regexpr2 = ...
    ['^\*(diff\(\w+(\(\W?\w+\))?\))' ...
    '|' '^\*\((\w+(\(\W?\w+\))?(\W?\w+(\(\W?\w+\))?)*)' ...
    '|' '^\*(\w+(\(\W?\w+\))?)'
    ];

M_endo_names_trim = cellfun(@strtrim, num2cell(M_.endo_names(:,:),2), 'Uniform', 0);
idxs = sortrows([(1:length(M_endo_names_trim))' cellfun(@length, M_endo_names_trim)], 2, 'descend');
regex = strjoin(M_endo_names_trim(idxs(:,1)), '|');
for i = 1:length(lhs)
    %% Construct regression matrices
    rhs_ = strsplit(rhs{i}, {'+','-','*','/','^','log(','exp(','(',')'});
    rhs_(cellfun(@(x) all(isstrprop(x, 'digit')), rhs_)) = [];
    vnames = setdiff(rhs_, cellstr(M_.param_names));
    if ~isempty(regexp(rhs{i}, ...
            ['(' strjoin(vnames, '\\(\\d+\\)|') '\\(\\d+\\))'], ...
            'once'))
        error(['dyn_ols: you cannot have leads in equation on line ' ...
            lineno{i} ': ' lhs{i} ' = ' rhs{i}]);
    end

    pnames = intersect(rhs_, cellstr(M_.param_names));
    vnames = cell(1, length(pnames));
    X = dseries();
    for j = 1:length(pnames)
        rhs_split = strsplit(rhs{i}, pnames{j});
        assert(length(rhs_split) == 2);
        if ~isempty(rhs_split{1}) && rhs_split{1}(end) == '*'
            tmp = regexp(rhs_split{1}, regexpr1, 'tokens');
        elseif ~isempty(rhs_split{2}) && rhs_split{2}(1) == '*'
            tmp = regexp(rhs_split{2}, regexpr2, 'tokens');
        else
            error('dyn_ols: Shouldn''t arrive here');
        end
        vnames{j} = tmp{1}{:};
        Xtmp = getdata(ds, regex, vnames{j});
        Xtmp.rename_(vnames{j});
        X = [X Xtmp];
    end
    Y = getdata(ds, regex, lhs{i}) ;

    fp = max(Y.firstobservedperiod, X.firstobservedperiod);
    lp = min(Y.lastobservedperiod, X.lastobservedperiod);

    Y = Y(fp:lp).data;
    X = X(fp:lp).data;

    %% Estimation
    % From LeSage, James P. "Applied Econometrics using MATLAB"
    if nargin == 2
        if iscell(varargin{1})
            tagv = varargin{1}{i};
        else
            tagv = varargin{1};
        end
    else
        tagv = ['eqlineno' num2str(lineno{i})];
    end
    [nobs, nvars] = size(X);
    oo_.ols.(tagv).dof = nobs - nvars;

    % Estimated Parameters
    [q, r] = qr(X, 0);
    xpxi = (r'*r)\eye(nvars);
    oo_.ols.(tagv).beta = r\(q'*Y);
    for j = 1:length(pnames)
        M_.params(strmatch(pnames{j}, M_.param_names, 'exact')) = oo_.ols.(tagv).beta(j);
    end

    % Yhat
    oo_.ols.(tagv).Yhat = X*oo_.ols.(tagv).beta;

    % Residuals
    oo_.ols.(tagv).resid = Y - oo_.ols.(tagv).Yhat;

    %% Calculate statistics
    % Estimate for sigma^2
    SS_res = oo_.ols.(tagv).resid'*oo_.ols.(tagv).resid;
    oo_.ols.(tagv).s2 = SS_res/oo_.ols.(tagv).dof;

    % R^2
    ym = Y - mean(Y);
    SS_tot = ym'*ym;
    oo_.ols.(tagv).R2 = 1 - SS_res/SS_tot;

    % Adjusted R^2
    oo_.ols.(tagv).adjR2 = oo_.ols.(tagv).R2 - (1 - oo_.ols.(tagv).R2)*nvars/(oo_.ols.(tagv).dof-1);

    % Durbin-Watson
    ediff = oo_.ols.(tagv).resid(2:nobs) - oo_.ols.(tagv).resid(1:nobs-1);
    oo_.ols.(tagv).dw = (ediff'*ediff)/SS_res;

    % Standard Error
    oo_.ols.(tagv).stderr = sqrt(oo_.ols.(tagv).s2*diag(xpxi));

    % T-Stat
    oo_.ols.(tagv).tstat = oo_.ols.(tagv).beta./oo_.ols.(tagv).stderr;

    %% Print Output
    title = sprintf('OLS Estimation of equation  `%s`', tagv);
    if nargin == 3
        title = [title sprintf(' [%s = %s]', 'name', tagv)];
    end

    preamble = {sprintf('Dependent Variable: %s', lhs{i}), ...
        sprintf('No. Independent Variables: %d', nvars), ...
        sprintf('Observations: %d', nobs)};

    afterward = {sprintf('R^2: %f', oo_.ols.(tagv).R2), ...
        sprintf('R^2 Adjusted: %f', oo_.ols.(tagv).adjR2), ...
        sprintf('s^2: %f', oo_.ols.(tagv).s2), ...
        sprintf('Durbin-Watson: %f', oo_.ols.(tagv).dw)};

    dyn_table(title, preamble, afterward, vnames, ...
        {'Coefficients','t-statistic','Std. Error'}, 4, ...
        [oo_.ols.(tagv).beta oo_.ols.(tagv).tstat oo_.ols.(tagv).stderr]);
end
end

function retval = getdata(ds, regex, ser)
if strncmp(ser, 'diff', 4)
    ser = ser(6:end-1);
    lagidx = strfind(ser, '(');
    if isempty(lagidx)
        retval = ds{ser} - ds{ser}(-1);
    else
        lag = str2double(ser(lagidx+1:strfind(ser, ')')-1));
        assert(lag < 0);
        ser = ser(1:lagidx-1);
        retval = ds{ser}(lag) - ds{ser}(lag-1);
    end
else
    retval = eval(regexprep(ser, regex, 'ds.$&'));
end
end