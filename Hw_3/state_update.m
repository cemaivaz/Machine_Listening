function [ lup ] = state_update(obs, log_p)
%STATE_UPDATE Summary of this function goes here
%   Detailed explanation goes here

lup = log(obs(:)) + log_p;
end

