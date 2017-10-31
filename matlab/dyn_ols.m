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
M_endo_exo_names_trim = cellfun(@strtrim, ...
    [num2cell(M_.endo_names(:,:),2) ; num2cell(M_.exo_names(:,:),2)], ...
    'Uniform', 0);
regex = strjoin(M_endo_exo_names_trim(:,1), '|');
mathops = '[\+\*\^\-\/]';
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
        createdvar = false;
        pregex = [...
            mathops pnames{j} mathops ...
            '|^' pnames{j} mathops ...
            '|' mathops pnames{j} '$' ...
            ];
        [startidx, endidx] = regexp(rhs{i}, pregex, 'start', 'end');
        assert(length(startidx) == 1);
        if rhs{i}(startidx) == '*'
            vnames{j} = getStrMoveLeft(rhs{i}(1:startidx-1));
        elseif rhs{i}(endidx) == '*'
            vnames{j} = getStrMoveRight(rhs{i}(endidx+1:end));
        elseif rhs{i}(startidx) == '+' ...
                || rhs{i}(startidx) == '-' ...
                || rhs{i}(endidx) == '+' ...
                || rhs{i}(endidx) == '-'
            % intercept
            createdvar = true;
            if any(strcmp(M_endo_exo_names_trim, 'intercept'))
                [~, vnames{j}] = fileparts(tempname);
                vnames{j} = ['intercept_' vnames{j}];
                assert(~any(strcmp(M_endo_exo_names_trim, vnames{j})));
            else
                vnames{j} = 'intercept';
            end
        else
            error('dyn_ols: Shouldn''t arrive here');
        end
        if createdvar
            Xtmp = dseries(ones(ds.nobs, 1), ds.firstdate, vnames{j});
        else
            Xtmp = eval(regexprep(vnames{j}, regex, 'ds.$&'));
            Xtmp.rename_(vnames{j});
        end
        X = [X Xtmp];
    end
    Y = eval(regexprep(lhs{i}, regex, 'ds.$&'));

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

function retval = getStrMoveLeft(str)
mathops = '[\+\*\^\-\/]';
if str(end) ~= ')'
    retval = str(max(regexp(str, mathops))+1:end);
else
    closedidxs = strfind(str, ')');
    closedidxs = [(length(closedidxs):-1:1)' closedidxs'];
    openidxs = strfind(str, '(');
    openidxs = [(length(openidxs):-1:1)' openidxs'];
    assert(rows(closedidxs) == rows(openidxs));
    for i = rows(openidxs):-1:1
        openparenidx = find(openidxs(i, 2) < closedidxs(:, 2), 1, 'first');
        if openidxs(i, 1) == closedidxs(openparenidx, 1)
            break;
        end
    end
    retval = str(max(regexp(str(1:openidxs(openparenidx,2)), mathops))+1:closedidxs(end));
end
end

function retval = getStrMoveRight(str)
mathops = '[\+\*\^\-\/]';
mathidxs = regexp(str, mathops);
openidxs = strfind(str, '(');
openidxs = [(1:length(openidxs))' openidxs'];
if min(mathidxs) < min(openidxs(:, 2))
    retval = str(1:min(regexp(str, mathops))-1);
else
    closedidxs = strfind(str, ')');
    closedidxs = [(1:length(closedidxs))' closedidxs'];
    assert(length(openidxs) == length(closedidxs));
    for i = 1:length(closedidxs)
        closedparenidx = sum(openidxs(:, 2) < closedidxs(i, 2));
        if openidxs(closedparenidx, 1) == closedidxs(i, 1)
            break;
        end
    end
    retval = str(1:closedidxs(closedparenidx, 2));
end
end
