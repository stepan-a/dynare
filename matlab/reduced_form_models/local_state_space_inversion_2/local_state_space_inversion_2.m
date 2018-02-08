function epsilon = local_state_space_inversion_2(y, yhat, ghx, ghu, constant, ghxx, ghuu, ghxu, numthread)

% Given the states (y) and structural innovations (epsilon), this routine computes the level of selected endogenous variables when the
% model is approximated by an order two taylor expansion around the deterministic steady state. Depending on the number of input/output
% argument the pruning algorithm advocated by C. Sims is used or not (this version should not be used if the selected endogenous variables
% are not the states of the model).
%
% INPUTS 
% - y          [double]           q*1 vector,
% - yhat       [double]           n*1 vector, initial conditions (n is the number of state variables).
% - ghx        [double]           q*n matrix, restricted dr.ghx where we only the lines corresponding to a subset of endogenous variables.
% - ghu        [double]           q*q matrix, restricted dr.ghu where we only consider the lines corresponding to a subset of endogenous variables.
% - constant   [double]           q*1 vector, deterministic steady state plus second order correction for a subset of endogenous variables.
% - ghxx       [double]           q*n² matrix, restricted dr.ghxx where we only consider the lines corresponding to a subset of endogenous variables.
% - ghuu       [double]           q*q² matrix, restricted dr.ghuu where we only consider the lines corresponding to a subset of endogenous variables.
% - ghxu       [double]           q*(nq) matrix, restricted dr.ghxu where we only consider the lines corresponding to a subset of endogenous variables.
% - numthead   [integer]          scalar, number of threads if parallelized dlls are available.
%
% OUTPUTS 
% - epsilon    [double]           q*1 vector, structural innovations.
%
% REMARKS 
% - [1] We do not consider the case with pruning. 
% - [2] 

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

q = rows(y);
n = rows(yhat);
s = columns(yhat);

if ~isequal(rows(ghx), q) || ~isequal(columns(ghx), n)
    error('Inconsitent dimensions between y, yhat and ghx')
end

if ~isequal(rows(ghu), q) || ~isequal(columns(ghu), q)
    error('Inconsitent dimensions between y and ghu')
end

if ~isequal(rows(ghxx), q) || ~isequal(columns(ghxx), n*n)
    error('Inconsitent dimensions between y, yhat and ghxx')
end

if ~isequal(rows(ghuu), q) || ~isequal(columns(ghuu), q*q)
    error('Inconsitent dimensions between y, and ghuu')
end

if ~isequal(rows(ghxu), q) || ~isequal(columns(ghxu), n*q)
    error('Inconsitent dimensions between y, yhat and ghxu')
end

if ~isequal(columns(y), s)
    error('y and yhat must have the same number of columns!')
end

epsilon = zeros(q, s);

for t = 1:s    
    epsilon0 = ghu\(y(:,t) - constant - ghx*yhat(:,t) - A_times_B_kronecker_C(.5*ghxx,yhat(:,t), numthread));
    iteration = 1;
    while iteration<=100
        f = y(:,t) - constant - ghx*yhat(:,t) - ghu*epsilon0 ...
            - A_times_B_kronecker_C(.5*ghxx, yhat(:,t), numthread)  ...
            - A_times_B_kronecker_C(.5*ghuu, epsilon0, numthread) ...
            - A_times_B_kronecker_C(ghxu, yhat(:,t), epsilon0, numthread);
        if f'*f<1e-5
            epsilon(:,t) = epsilon0;
            break
        end
        df = -(ghu + .5*ghuu*(kron(eye(q), epsilon0)+kron(epsilon0, eye(q))) + ghxu*kron(yhat(:,t), eye(q)));
        epsilon1 = epsilon0 - df\f;
        if (epsilon1-epsilon0)'*(epsilon1-epsilon0)<1e-5
            epsilon(:,t) = epsilon1;
            break
        end
        iteration = iteration+1;
    end
    if ~all(epsilon(:,t))
        error('Newton did not converge in period %s', num2str(t))
    end
end