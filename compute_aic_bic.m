function [AIC,BIC,LL] = compute_aic_bic (qc,model,params)
% Compute AIC and BIC for bespoke multinomial model
% FORMAT [AIC,BIC,LL] = compute_aic_bic (qc,model,params)
%
% qc        behaviour correct rates
% model     model correct rates
% params    numbers of parameters and trials
%
% LL        Log Likelihood

names=params.names;
T=params.T;
k=params.k;

% Number of data points
N = sum(T);

q=qc/100;

p_delta=0.01;
%p_delta=eps;

for m=1:length(model),
    model(m).p=model(m).pc/100; % Probabilities
    
    % Log Likelihoods (under binomial model)
    like(m,:)=q.*log(model(m).p+p_delta) + (1-q).*log(1-model(m).p+p_delta);
    LL(m) = sum(T.*like(m,:));
    
    e=q-model(m).p;
    MSE(m)=mean(e.^2);
end

debug=0;
if debug
    disp(like);
    keyboard
end

AIC = -2*LL+2*k;
BIC = -2*LL+k*log(N);

disp('AIC and BIC : less positive the better');
figure
subplot(2,2,1);
bar(LL);
grid on
set(gca,'XTickLabel',names)
ylabel('Log Likelihood');
subplot(2,2,2);
bar(AIC);
grid on
set(gca,'XTickLabel',names)
ylabel('AIC');
subplot(2,2,4);
bar(BIC);
grid on
set(gca,'XTickLabel',names)
ylabel('BIC');
subplot(2,2,3);
bar(MSE);
grid on
set(gca,'XTickLabel',names)
ylabel('MSE');

