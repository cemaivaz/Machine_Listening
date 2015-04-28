function [ lpp ] = state_predict(A, log_p)
%LOG_ALPHA_PREDICT Summary of this function goes here
%   Detailed explanation goes here
mx = max(log_p(:));
p = exp(log_p - mx);
lpp = log(A * p) + mx;

end

