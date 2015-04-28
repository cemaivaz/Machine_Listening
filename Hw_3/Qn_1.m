%% setup parameter structure

K = 200; %length of sequence, e.g. at sampling period of 100ms: 50s 
N = 16; %number of states: how many dicisions has a bar?
eps = 0.6; %probability to advance, higher value: higher tempo
nu = 0.9; %probability to observe onset on downbeat
%flat prior 
pri = ones(N,1)/N;
%transition matrix: self-transition, or one step ahead.
I = eye(N);
A = (1-eps)*I + eps*I(:, [2:N, 1]);
%observation model: very likely to observe onsets on the downbeat,
% and very unlikely elsewhere: very simplistic dummy model...
C = [nu nu/2 (1-nu)*ones(1,N-2)];
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

beta = zeros(N, K);
beta_postdict = zeros(N, K);
for t=K:-1:1,
    if t==K,
        beta_postdict(:,t) = ones(N,1);
    else
        beta_postdict(:,t) = hm.A'*beta(:,t+1);
    end;
    beta(:, t) = hm.C(obs(t)+1, :)' .* beta_postdict(:,t);
    beta(:,t) = beta(:,t)./sum(beta(:,t));
end;

subplot(211)
imagesc(beta)
hold on
plot(state,'w--','LineWidth',2);
title('Update messages','FontSize',20)
set(gca,'FontSize',20)
subplot(212)
imagesc(beta_postdict)
hold on
plot(state,'w--','LineWidth',2);
title('Postdict messages','FontSize',20)
set(gca,'FontSize',20)

%% Smoothing

gamma = log_alpha .* beta_postdict; 

subplot(211)
hold off
imagesc(gamma)
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