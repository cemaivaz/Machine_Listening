K = 200; %length of sequence, e.g. at sampling period of 100ms: 50s
N = 10; %number of states: how many divisions has a bar?
eps = 0.6; %probability to advance, higher value: higher tempo
nu = 0.9; %probability to observe onset on downbeat
%flat prior
pri = ones(N,1)/N;
%transition matrix: self-transition, or one step ahead.
I = eye(N);
A = (1-eps)*I + eps*I(:, [2:N, 1]);
%observation model: assumes high likelihood for observing...
%...an onset at a stroke of the Curcuna usul.
C = [nu 1-nu nu nu 1-nu nu 1-nu nu 1-nu nu];%[nu nu/2 (1-nu)*ones(1,N-2)];
C = [1-C;C];
%data structure for the HMM model!
hm = struct('A', A, 'C', C, 'p_x1', pri);

state = zeros(1, K);
obs = zeros(1, K);
for t=1:K,
if t==1,
state(t)=randsample(1:N,1,true, hm.p_x1);
else
state(t)=randsample(1:N,1,true, hm.A(:, state(t-1)));
end
obs(t) = randsample([0,1], 1, true, C(:, state(t)));
end
%%

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
%%

alpha_pred = zeros(N, K);
alpha = zeros(N, K);
for t=1:K,
if t==1,
alpha_pred(:,t) = pri;
else
alpha_pred(:,t) = A*alpha(:,t-1);
end
alpha(:,t) = alpha_pred(:,t).*(C(obs(t)+1,:)');
alpha(:,t) = alpha(:,t)./sum(alpha(:,t));
%note: normalization enables visualization. Alphas
%and betas decrease throughout recursions.
end

subplot(211)
imagesc(alpha)
hold on
plot(state,'w--','LineWidth',2);
title('Update messages','FontSize',20)
set(gca,'FontSize',20)
subplot(212)
imagesc(alpha_pred)
hold on
plot(state,'w--','LineWidth',2);
title('Predict messages','FontSize',20)
set(gca,'FontSize',20)

%%
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
%%
gamma = alpha .* beta_postdict;



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