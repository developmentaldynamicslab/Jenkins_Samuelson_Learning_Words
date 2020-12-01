%Plots mean results across subjects for Bayes simulations.  Input current
%peakHistory array followed by the arrays containing positions of all items
%in all sets.

function [graphMatrix] = plotBayes( peakHistory, trialtype, simultaneous) 

%behavioral data to compare to (condition x test type[sub bas sup])
E1avg(1,1) = 99.12263158;
E1avg(2,1) = 98.24526316;
E1avg(3,1) = 96.49052632;
E1avg(4,1) = 94.73631579;
E1avg(1,2) = 48.24263158;
E1avg(2,2) = 10.52526316;
E1avg(3,2) = 92.10421053;
E1avg(4,2) = 85.08631579;
E1avg(1,3) = 7.016315789;
E1avg(2,3) = 0.876842105;
E1avg(3,3) = 15.34684211;
E1avg(4,3) = 81.31210526;
E2avg(1,1) = 94.9995;
E2avg(2,1) = 88.3315;
E2avg(3,1) = 94.1665;
E2avg(4,1) = 90.833;
E2avg(1,2) = 30.832;
E2avg(2,2) = 53.3305;
E2avg(3,2) = 89.9985;
E2avg(4,2) = 86.6655;
E2avg(1,3) = 5.8325;
E2avg(2,3) = 2.499;
E2avg(3,3) = 13.33;
E2avg(4,3) = 75.412;
E3avg(1,1) = 94.73631579;
E3avg(2,1) = 92.98210526;
E3avg(3,1) = 98.24526316;
E3avg(4,1) = 93.85894737;
E3avg(1,2) = 39.90947368;
E3avg(2,2) = 51.75315789;
E3avg(3,2) = 92.10421053;
E3avg(4,2) = 72.80578947;
E3avg(1,3) = 4.824210526;
E3avg(2,3) = 14.47263158;
E3avg(3,3) = 17.54105263;
E3avg(4,3) = 67.97947368;

sH = size(peakHistory);
sH = sH(2) / 2;
subjects = max(peakHistory(9,:));
rounds = max(peakHistory(10,:));
trials = max(peakHistory(11,:));

sngPercSub = 0;
sngPercBas = 0;
sngPercSup = 0;
tsubPercSub = 0;
tsubPercBas = 0;
tsubPercSup = 0;
tbasPercSub = 0;
tbasPercBas = 0;
tbasPercSup = 0;
tsupPercSub = 0;
tsupPercBas = 0;
tsupPercSup = 0;
singleRuns = 0;
tsubRuns = 0;
tbasRuns = 0;
tsupRuns = 0;

sngPercSubS = [];
sngPercBasS = [];
sngPercSupS = [];
tsubPercSubS = [];
tsubPercBasS = [];
tsubPercSupS = [];
tbasPercSubS = [];
tbasPercBasS = [];
tbasPercSupS = [];
tsupPercSubS = [];
tsupPercBasS = [];
tsupPercSupS = [];

for ss = 1:subjects
    for rs = 1:rounds
        for ts = 1:trials
            currsims = ts*2 + (rs-1)*trials*2 + (ss-1)*rounds*trials*2 - 1; %this brings you to the columns that say 1 or 0 based on if a peak was there.
            setType = mod(rs,3) + 1;%animals = 1, vehicles = 2, vegetables = 3
            if peakHistory(12,currsims) == 1 %trial type = single
                singleRuns = singleRuns + 1; %used for total possible
                if trialtype == 4
                    sngPercBasS = [sngPercBasS peakHistory(1,currsims)];
                    sngPercBasS = [sngPercBasS peakHistory(2,currsims)];
                    if peakHistory(1,currsims) == 1
                        sngPercBas = sngPercBas + 1;
                    end
                    if peakHistory(2,currsims) == 1
                        sngPercBas = sngPercBas + 1;
                    end
                elseif trialtype == 5
                    sngPercSubS = [sngPercSubS peakHistory(1,currsims)];
                    sngPercSubS = [sngPercSubS peakHistory(2,currsims)];
                    if peakHistory(1,currsims) == 1
                        sngPercSub = sngPercSub + 1;
                    end
                    if peakHistory(2,currsims) == 1
                        sngPercSub = sngPercSub + 1;
                    end
                elseif trialtype == 6
                    sngPercBasS = [sngPercBasS peakHistory(1,currsims)];
                    sngPercBasS = [sngPercBasS peakHistory(2,currsims)];
                    if peakHistory(1,currsims) == 1
                        sngPercBas = sngPercBas + 1;
                    end
                    if peakHistory(2,currsims) == 1
                        sngPercBas = sngPercBas + 1;
                    end
                elseif trialtype == 7
                    sngPercSupS = [sngPercSupS peakHistory(1,currsims)];
                    sngPercSupS = [sngPercSupS peakHistory(2,currsims)];
                    sngPercSupS = [sngPercSupS peakHistory(3,currsims)];
                    sngPercSupS = [sngPercSupS peakHistory(4,currsims)];
                    if peakHistory(1,currsims) == 1
                        sngPercSup = sngPercSup + 1;
                    end
                    if peakHistory(2,currsims) == 1
                        sngPercSup = sngPercSup + 1;
                    end
                    if peakHistory(3,currsims) == 1
                        sngPercSup = sngPercSup + 1;
                    end
                    if peakHistory(4,currsims) == 1
                        sngPercSup = sngPercSup + 1;
                    end
                else
                    for sing = 1:8 %used for total observed
                        if (sing == 1 || sing == 5)
                            sngPercSubS = [sngPercSubS peakHistory(sing,currsims)];
                        elseif  (sing == 3 || sing == 8)
                            sngPercBasS = [sngPercBasS peakHistory(sing,currsims)];
                        else
                            sngPercSupS = [sngPercSupS peakHistory(sing,currsims)];
                        end
                        
                        if (sing == 1 || sing == 5) && peakHistory(sing,currsims) == 1
                            sngPercSub = sngPercSub + 1;
                        elseif  (sing == 3 || sing == 8) && peakHistory(sing,currsims) == 1
                            sngPercBas = sngPercBas + 1;
                        elseif  peakHistory(sing,currsims) == 1
                            sngPercSup = sngPercSup + 1;
                        end
                    end
                end
            end
            if peakHistory(12,currsims) == 2 %trial type = tsub
                tsubRuns = tsubRuns + 1; %used for total possible
                if trialtype == 4
                    tsubPercBasS = [tsubPercBasS peakHistory(1,currsims)];
                    tsubPercBasS = [tsubPercBasS peakHistory(2,currsims)];
                    if peakHistory(1,currsims) == 1
                        tsubPercBas = tsubPercBas + 1;
                    end
                    if peakHistory(2,currsims) == 1
                        tsubPercBas = tsubPercBas + 1;
                    end
                elseif trialtype == 5
                    tsubPercSubS = [tsubPercSubS peakHistory(1,currsims)];
                    tsubPercSubS = [tsubPercSubS peakHistory(2,currsims)];
                    if peakHistory(1,currsims) == 1
                        tsubPercSub = tsubPercSub + 1;
                    end
                    if peakHistory(2,currsims) == 1
                        tsubPercSub = tsubPercSub + 1;
                    end
                elseif trialtype == 6
                    tsubPercBasS = [tsubPercBasS peakHistory(1,currsims)];
                    tsubPercBasS = [tsubPercBasS peakHistory(2,currsims)];
                    if peakHistory(1,currsims) == 1
                        tsubPercBas = tsubPercBas + 1;
                    end
                    if peakHistory(2,currsims) == 1
                        tsubPercBas = tsubPercBas + 1;
                    end
                elseif trialtype == 7
                    tsubPercSupS = [tsubPercSupS peakHistory(1,currsims)];
                    tsubPercSupS = [tsubPercSupS peakHistory(2,currsims)];
                    tsubPercSupS = [tsubPercSupS peakHistory(3,currsims)];
                    tsubPercSupS = [tsubPercSupS peakHistory(4,currsims)];
                    if peakHistory(1,currsims) == 1
                        tsubPercSup = tsubPercSup + 1;
                    end
                    if peakHistory(2,currsims) == 1
                        tsubPercSup = tsubPercSup + 1;
                    end
                    if peakHistory(3,currsims) == 1
                        tsubPercSup = tsubPercSup + 1;
                    end
                    if peakHistory(4,currsims) == 1
                        tsubPercSup = tsubPercSup + 1;
                    end
                else
                    for tsub = 1:8 %used for total observed
                        if (tsub == 1 || tsub == 5)
                            tsubPercSubS = [tsubPercSubS peakHistory(tsub,currsims)];
                        elseif  (tsub == 3 || tsub == 8)
                            tsubPercBasS = [tsubPercBasS peakHistory(tsub,currsims)];
                        else
                            tsubPercSupS = [tsubPercSupS peakHistory(tsub,currsims)];
                        end
                        
                        if (tsub == 1 || tsub == 5) && peakHistory(tsub,currsims) == 1
                            tsubPercSub = tsubPercSub + 1;
                        elseif (tsub == 3 || tsub == 8) && peakHistory(tsub,currsims) == 1
                            tsubPercBas = tsubPercBas + 1;
                        elseif  peakHistory(tsub,currsims) == 1
                            tsubPercSup = tsubPercSup + 1;
                        end
                    end
                end
            end
            if peakHistory(12,currsims) == 3 %trial type = tbas
                tbasRuns = tbasRuns + 1; %used for total possible
                if trialtype == 4
                    tbasPercBasS = [tbasPercBasS peakHistory(1,currsims)];
                    tbasPercBasS = [tbasPercBasS peakHistory(2,currsims)];
                    if peakHistory(1,currsims) == 1
                        tbasPercBas = tbasPercBas + 1;
                    end
                    if peakHistory(2,currsims) == 1
                        tbasPercBas = tbasPercBas + 1;
                    end
                elseif trialtype == 5
                    tbasPercSubS = [tbasPercSubS peakHistory(1,currsims)];
                    tbasPercSubS = [tbasPercSubS peakHistory(2,currsims)];
                    if peakHistory(1,currsims) == 1
                        tbasPercSub = tbasPercSub + 1;
                    end
                    if peakHistory(2,currsims) == 1
                        tbasPercSub = tbasPercSub + 1;
                    end
                elseif trialtype == 6
                    tbasPercBasS = [tbasPercBasS peakHistory(1,currsims)];
                    tbasPercBasS = [tbasPercBasS peakHistory(2,currsims)];
                    if peakHistory(1,currsims) == 1
                        tbasPercBas = tbasPercBas + 1;
                    end
                    if peakHistory(2,currsims) == 1
                        tbasPercBas = tbasPercBas + 1;
                    end
                elseif trialtype == 7
                    tbasPercSupS = [tbasPercSupS peakHistory(1,currsims)];
                    tbasPercSupS = [tbasPercSupS peakHistory(2,currsims)];
                    tbasPercSupS = [tbasPercSupS peakHistory(3,currsims)];
                    tbasPercSupS = [tbasPercSupS peakHistory(4,currsims)];
                    if peakHistory(1,currsims) == 1
                        tbasPercSup = tbasPercSup + 1;
                    end
                    if peakHistory(2,currsims) == 1
                        tbasPercSup = tbasPercSup + 1;
                    end
                    if peakHistory(3,currsims) == 1
                        tbasPercSup = tbasPercSup + 1;
                    end
                    if peakHistory(4,currsims) == 1
                        tbasPercSup = tbasPercSup + 1;
                    end
                else
                    for tbas = 1:8 %used for total observed
                        if (tbas == 1 || tbas == 5)
                            tbasPercSubS = [tbasPercSubS peakHistory(tbas,currsims)];
                        elseif  (tbas == 3 || tbas == 8)
                            tbasPercBasS = [tbasPercBasS peakHistory(tbas,currsims)];
                        else
                            tbasPercSupS = [tbasPercSupS peakHistory(tbas,currsims)];
                        end
                        
                        if (tbas == 1 || tbas == 5) && peakHistory(tbas,currsims) == 1
                            tbasPercSub = tbasPercSub + 1;
                        elseif (tbas == 3 || tbas == 8) && peakHistory(tbas,currsims) == 1
                            tbasPercBas = tbasPercBas + 1;
                        elseif peakHistory(tbas,currsims) == 1
                            tbasPercSup = tbasPercSup + 1;
                        end
                    end
                end
            end
            if peakHistory(12,currsims) == 4 %trial type = tsup
                tsupRuns = tsupRuns + 1; %used for total possible
                if trialtype == 4
                    tsupPercBasS = [tsupPercBasS peakHistory(1,currsims)];
                    tsupPercBasS = [tsupPercBasS peakHistory(2,currsims)];
                    if peakHistory(1,currsims) == 1
                        tsupPercBas = tsupPercBas + 1;
                    end
                    if peakHistory(2,currsims) == 1
                        tsupPercBas = tsupPercBas + 1;
                    end
                elseif trialtype == 5
                    tsupPercSubS = [tsupPercSubS peakHistory(1,currsims)];
                    tsupPercSubS = [tsupPercSubS peakHistory(2,currsims)];
                    if peakHistory(1,currsims) == 1
                        tsupPercSub = tsupPercSub + 1;
                    end
                    if peakHistory(2,currsims) == 1
                        tsupPercSub = tsupPercSub + 1;
                    end
                elseif trialtype == 6
                    tsupPercBasS = [tsupPercBasS peakHistory(1,currsims)];
                    tsupPercBasS = [tsupPercBasS peakHistory(2,currsims)];
                    if peakHistory(1,currsims) == 1
                        tsupPercBas = tsupPercBas + 1;
                    end
                    if peakHistory(2,currsims) == 1
                        tsupPercBas = tsupPercBas + 1;
                    end
                elseif trialtype == 7
                    tsupPercSupS = [tsupPercSupS peakHistory(1,currsims)];
                    tsupPercSupS = [tsupPercSupS peakHistory(2,currsims)];
                    tsupPercSupS = [tsupPercSupS peakHistory(3,currsims)];
                    tsupPercSupS = [tsupPercSupS peakHistory(4,currsims)];
                    if peakHistory(1,currsims) == 1
                        tsupPercSup = tsupPercSup + 1;
                    end
                    if peakHistory(2,currsims) == 1
                        tsupPercSup = tsupPercSup + 1;
                    end
                    if peakHistory(3,currsims) == 1
                        tsupPercSup = tsupPercSup + 1;
                    end
                    if peakHistory(4,currsims) == 1
                        tsupPercSup = tsupPercSup + 1;
                    end
                else
                    for tsup = 1:8 %used for total observed
                        if (tsup == 1 || tsup == 5)
                            tsupPercSubS = [tsupPercSubS peakHistory(tsup,currsims)];
                        elseif  (tsup == 3 || tsup == 8)
                            tsupPercBasS = [tsupPercBasS peakHistory(tsup,currsims)];
                        else
                            tsupPercSupS = [tsupPercSupS peakHistory(tsup,currsims)];
                        end
                        
                        if (tsup == 1 || tsup == 5) && peakHistory(tsup,currsims) == 1
                            tsupPercSub = tsupPercSub + 1;
                        elseif (tsup == 3 || tsup == 8) && peakHistory(tsup,currsims) == 1
                            tsupPercBas = tsupPercBas + 1;
                        elseif peakHistory(tsup,currsims) == 1
                            tsupPercSup = tsupPercSup + 1;
                        end
                    end
                end
            end
        end
    end
end

%Convert to percentages
if singleRuns > 0
    sngPercSub = sngPercSub / (singleRuns*2);
    sngPercBas = sngPercBas / (singleRuns*2);
    sngPercSup = sngPercSup / (singleRuns*4);
end
if tsubRuns > 0
    tsubPercSub = tsubPercSub / (tsubRuns*2);
    tsubPercBas = tsubPercBas / (tsubRuns*2);
    tsubPercSup = tsubPercSup / (tsubRuns*4);
end
if tbasRuns > 0
    tbasPercSub = tbasPercSub / (tbasRuns*2);
    tbasPercBas = tbasPercBas / (tbasRuns*2);
    tbasPercSup = tbasPercSup / (tbasRuns*4);
end
if tsupRuns > 0
    tsupPercSub = tsupPercSub / (tsupRuns*2);
    tsupPercBas = tsupPercBas / (tsupRuns*2);
    tsupPercSup = tsupPercSup / (tsupRuns*4);
end


%plot
graphMatrix = [sngPercSub sngPercBas sngPercSup;tsubPercSub tsubPercBas tsubPercSup;tbasPercSub tbasPercBas tbasPercSup;tsupPercSub tsupPercBas tsupPercSup];
figure;
bar(graphMatrix);
axis([0 5 0 1]);
xlabel('   Sng     Sub    Bas     Sup');
ylabel('% Chosen');
legend('sub','bas','sup');
if(simultaneous == 1)
    title('simultaneous');
elseif(simultaneous == 0)
    title('sequential');
elseif(simultaneous == 2)
    title('experiment 3');
end

disp('plot should appear.');
sngPercSubS = sqrt(((1/(size(sngPercSubS,2)-1)) * sum((sngPercSubS - mean(sngPercSubS)).^2))) / sqrt(size(sngPercSubS,2));
sngPercBasS = sqrt(((1/(size(sngPercBasS,2)-1)) * sum((sngPercBasS - mean(sngPercBasS)).^2))) / sqrt(size(sngPercBasS,2));
sngPercSupS = sqrt(((1/(size(sngPercSupS,2)-1)) * sum((sngPercSupS - mean(sngPercSupS)).^2))) / sqrt(size(sngPercSupS,2));
tsubPercSubS = sqrt(((1/(size(tsubPercSubS,2)-1)) * sum((tsubPercSubS - mean(tsubPercSubS)).^2))) / sqrt(size(tsubPercSubS,2));
tsubPercBasS = sqrt(((1/(size(tsubPercBasS,2)-1)) * sum((tsubPercBasS - mean(tsubPercBasS)).^2))) / sqrt(size(tsubPercBasS,2));
tsubPercSupS = sqrt(((1/(size(tsubPercSupS,2)-1)) * sum((tsubPercSupS - mean(tsubPercSupS)).^2))) / sqrt(size(tsubPercSupS,2));
tbasPercSubS = sqrt(((1/(size(tbasPercSubS,2)-1)) * sum((tbasPercSubS - mean(tbasPercSubS)).^2))) / sqrt(size(tbasPercSubS,2));
tbasPercBasS = sqrt(((1/(size(tbasPercBasS,2)-1)) * sum((tbasPercBasS - mean(tbasPercBasS)).^2))) / sqrt(size(tbasPercBasS,2));
tbasPercSupS = sqrt(((1/(size(tbasPercSupS,2)-1)) * sum((tbasPercSupS - mean(tbasPercSupS)).^2))) / sqrt(size(tbasPercSupS,2));
tsupPercSubS = sqrt(((1/(size(tsupPercSubS,2)-1)) * sum((tsupPercSubS - mean(tsupPercSubS)).^2))) / sqrt(size(tsupPercSubS,2));
tsupPercBasS = sqrt(((1/(size(tsupPercBasS,2)-1)) * sum((tsupPercBasS - mean(tsupPercBasS)).^2))) / sqrt(size(tsupPercBasS,2));
tsupPercSupS = sqrt(((1/(size(tsupPercSupS,2)-1)) * sum((tsupPercSupS - mean(tsupPercSupS)).^2))) / sqrt(size(tsupPercSupS,2));

OutName = ['Bayes3Result_' num2str(simultaneous) datestr(now, '_yyyy-mm-dd-THHMMSS') '.mat'];
save(OutName,'sngPercSub','sngPercBas','sngPercSup','tsubPercSub','tsubPercBas','tsubPercSup', ...
    'tbasPercSub','tbasPercBas','tbasPercSup','tsupPercSub','tsupPercBas','tsupPercSup', ...
    'sngPercSubS','sngPercBasS','sngPercSupS','tsubPercSubS','tsubPercBasS','tsubPercSupS', ...
    'tbasPercSubS','tbasPercBasS','tbasPercSupS','tsupPercSubS','tsupPercBasS','tsupPercSupS');

    
OutName2 = ['Bayes3Result_' num2str(simultaneous) datestr(now, '_yyyy-mm-dd-THHMMSS') '.prn'];
OutFile2 = fopen(OutName2,'a');
fprintf(OutFile2,'Exp N Sb1Sub Bs1Sub Sp1Sub Sb3Sub Bs3Sub Sp3Sub Sb3Bas Bs3Bas Sp3Bas Sb3Sup Bs3Sup Sp3Sup AIC-DFT AIC-Bayes AIC-Dummy BIC-DFT BIC-Bayes BIC-Dummy LL-DFT LL-Bayes LL-Dummy\n');

%New 2020 AIC/BIC
model(1).pc = 100.*[sngPercSub,sngPercBas,sngPercSup,tsubPercSub,tsubPercBas,tsubPercSup,tbasPercSub,tbasPercBas,tbasPercSup,tsupPercSub,tsupPercBas,tsupPercSup]; % DNF
if simultaneous == 1
    data = E1avg;
    model(2).pc = [98,78,12,89,29,1,100,100,8,100,100,100]; % Bayes
    model(3).pc = 50*ones(1,12); % Dummy
elseif simultaneous == 0
    data = E2avg;
    model(2).pc = [98,78,12,89,29,1,100,100,8,100,100,100]; % Bayes
    model(3).pc = 50*ones(1,12); % Dummy
elseif simultaneous == 2
    data = E3avg;
    model(2).pc = [98,78,12,86,2,0,100,100,3,100,100,100]; % Bayes
    model(3).pc = 50*ones(1,12); % Dummy
end

ct=1;
qc=zeros(1,12);
for i = 1:4 %1sub 3sub 3bas 3sup
    for k=1:3 %sub bas super
        qc(1,ct) = data(i,k);
        ct=ct+1;
    end
end

params.names={'DNF','Bayes','Dummy'};
% Number of model parameters
params.k=[19,2,0];
% Numbers of trials of each type
T = [6,6,12,6,6,12,6,6,12,6,6,12]; % per subject
params.T = T * 19; % all subjects

[AIC,BIC,LL] = compute_aic_bic (qc,model,params);

fprintf(OutFile2,'%i %i %6.2f %6.2f %6.2f %6.2f %6.2f %6.2f %6.2f %6.2f %6.2f %6.2f %6.2f %6.2f %6.2f %6.2f %6.2f %6.2f %6.2f %6.2f %6.2f %6.2f %6.2f\n', ...
    simultaneous,subjects,model(1).pc,AIC,BIC,LL);
    
fclose(OutFile2);


end

