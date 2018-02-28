function alpha = a2alpha(a)

% Computes the m alpha coefficients from the m a coefficients of the PAC model.
%
% INPUTS 
% - a      [double]   m*1 vector of coefficients.
%
% OUTPUTS 
% - alpha  [double]   m*1 vector of coefficients.
%
% NOTES 
%
%  Given the current estimate of the PAC parameters a_0, a_1, ..., a_{m-1}, the routine does the following:
%
%  \alpha_{m} = a_{m-1}
%  \alpha_{m-1} = a_{m-2}-a_{m-1}
%  \alpha_{m-2} = a_{m-3}-a_{m-2}
%  ...
%  \alpha_3 = a_2-a_3
%  \alpha_2 = a_1-a_2
%  \alpha_1 = a_0-a_1-1
%
% Note that the last elements of input a are (a_0, a_1, ..., a_{m-1}).

% Return an error if the input is not a vector
if ~isvector(a)
    error('Input argument has to be a vector of doubles!')
end 

% Get the number of PAC parameters (without the discount factor)
m = length(a);

% Initialize the vector of transformed PAC parameters.
alpha = zeros(m, 1);

% Compute the transformed parameters
alpha(m) = a(m);
alpha(2:m-1) = a(2:m-1)-a(3:m);
alpha(1) = a(1)-a(2)-1;