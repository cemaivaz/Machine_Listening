%% setup parameter structure

K = 200; %length of sequence, e.g. at sampling period of 100ms: 50s 
N = 10 * 2; %number of states: how many dicisions has a bar?
eps = 0.6; %probability to advance, higher value: higher tempo
nu = 0.9; %probability to observe onset on downbeat
%flat prior 
pri = ones(N,1)/N;
%transition matrix: self-transition, or one step ahead.
I = eye(N);
A = (1-eps)*I + eps*I(:, [2:N, 1]);
%observation model: very likely to observe onsets on the downbeat,
% and very unlikely elsewhere: very simplistic dummy model...
% C = [nu nu/2 (1-nu)*ones(1,N-2)];
% C = [1-C;C];
C = [0.9 0.088 0.27 0.097 0.72 0.086 0.81 0.092 0.4901 0.089 0.92 0.11 0.69 0.15 0.91 0.12 0.25 0.119 0.51 0.102];
C = [1-C;C];
%data structure for the HMM model!
hm = struct('A', A, 'C', C, 'p_x1', pri);

%% Generate data from the true model

state = zeros(1, K);
obs = zeros(1, K);
for t=1:K,
    if t==1,
        state(t) = randsample(1:N, 1, true, hm.p_x1);
    else
        state(t) = randsample(1:N, 1, true, hm.A(:, state(t-1)));
    end
    obs(t) = randsample([0,1], 1, true, C(:, state(t)));

end

%% plot the data:
figure;
subplot(211)
plot(state,'k');
title('Hidden state-sequence','FontSize',20)
set(gca,'FontSize',20)
subplot(212)
plot(obs,'xr-')
title('Observation sequence','FontSize',20)
set(gca,'FontSize',20)


%% forward pass

log_alpha_predict = zeros(N, K);
log_alpha = zeros(N, K);

for t=1:K,
    if t==1,
        log_alpha_predict(:,t) = log(hm.p_x1);
    else
%         log_alpha_predict(:,t) = A*log_alpha(:,t-1);
        %
        log_alpha_predict(:, t) = state_predict(hm.A, log_alpha(:, t-1));
    end
    log_alpha(:, t) = state_update(hm.C(obs(t) + 1, :), log_alpha_predict(:, t));
%     log_alpha(:,t) = log_alpha_predict(:,t).*(C(obs(t)+1,:)');
%     log_alpha(:,t) = log_alpha(:,t)./sum(log_alpha(:,t));


end

subplot(211)
imagesc(log_alpha)
hold on
plot(state,'w--','LineWidth',2);
title('Update messages','FontSize',20)
set(gca,'FontSize',20)
subplot(212)
imagesc(log_alpha_predict)
hold on
plot(state,'w--','LineWidth',2);
title('Predict messages','FontSize',20)
set(gca,'FontSize',20)

%% backward pass

log_beta = zeros(N, K);
log_beta_postdict = zeros(N, K);
for t=K:-1:1,
    if t==K,
        log_beta_postdict(:,t) = zeros(N,1);
        
    else
        log_beta_postdict(:, t) = state_postdict(hm.A, log_beta(:, t + 1));
%         log_beta_postdict(:,t) = hm.A'*log_beta(:,t+1);
    end;
    log_beta(:, t) = state_update(hm.C(obs(t) + 1, :), log_beta_postdict(:, t));
%     log_beta(:, t) = hm.C(obs(t)+1, :)' .* log_beta_postdict(:,t);
%     log_beta(:,t) = log_beta(:,t)./sum(log_beta(:,t));
end;

subplot(211)
imagesc(log_beta)
hold on
plot(state,'w--','LineWidth',2);
title('Update messages','FontSize',20)
set(gca,'FontSize',20)
subplot(212)
imagesc(log_beta_postdict)
hold on
plot(state,'w--','LineWidth',2);
title('Postdict messages','FontSize',20)
set(gca,'FontSize',20)

%% Smoothing

log_gamma = log_alpha + log_beta_postdict; 

subplot(211)
hold off
imagesc(log_gamma)
hold on
plot(state,'w--','LineWidth',2);
title('Inferred posterior','FontSize',20)
set(gca,'FontSize',20)
subplot(212)
hold off
plot(obs,'xr-')
title('Observations','FontSize',20)
set(gca,'FontSize',20)
colormap jet