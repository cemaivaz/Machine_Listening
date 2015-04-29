function [ lpp ] = state_postdict( A, log_p )
%STATE_POSTDICT Summary of this function goes here
%   Detailed explanation goes here
mx = max(log_p(:));
p = exp(log_p - mx);
lpp = log(A' * p) + mx;

end

