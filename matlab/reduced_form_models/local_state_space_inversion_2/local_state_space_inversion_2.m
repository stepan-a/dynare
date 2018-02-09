function [epsilon, exitflag] = local_state_space_inversion_2(y, yhat, ghx, ghu, constant, ghxx, ghuu, ghxu, numthread)

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

exitflag = 0;

epsilon  = zeros(q, 1);
epsilon0 = eguess(ghx, ghu, ghxx, constant, numthread, y, yhat);

iteration = 1;

while iteration<=100
    f0 = rf2r(ghx, ghu, ghxx, ghuu, ghxu, constant, numthread, y, yhat, epsilon0);
    if f0'*f0<1e-5
        epsilon = epsilon0;
        break
    end
    df0 = rf2J(ghu, ghuu, ghxu, yhat, epsilon0);
    epsilon1 = epsilon0 - df0\f0;
    f1 = rf2r(ghx, ghu, ghxx, ghuu, ghxu, constant, numthread, y, yhat, epsilon1);
    if (epsilon1-epsilon0)'*(epsilon1-epsilon0)<1e-5
        epsilon = epsilon1;
        break
    end
    epsilon0 = epsilon1;
    iteration = iteration+1;
end

if ~all(epsilon)
    exitflag = 1;
end