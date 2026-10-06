%% Fit the Tsinman-based PAMPA-BBB permeability model
% The model structure matches the constrained Collander descriptor used for
% the Smith et al. model:
%   log10(P0 [cm/s]) = intercept + slope*(logP - 0.5*log10(MW [g/mol]))

clear; clc;
scriptDir=fileparts(mfilename('fullpath'));
repoRoot=fileparts(fileparts(scriptDir));
inputFile=fullfile(repoRoot,'data','model_development','tsinman_2011_pampa_bbb.csv');
outputDir=fullfile(repoRoot,'results','model_development');
if ~isfile(inputFile), error('Tsinman dataset was not found: %s',inputFile); end
if ~exist(outputDir,'dir'), mkdir(outputDir); end

data=readtable(inputFile,'TextType','string');
required=["Compound","MW_g_per_mol","logP","log10_P0_cm_per_s"];
assert(all(ismember(required,string(data.Properties.VariableNames))), ...
    'The Tsinman dataset does not contain the required columns.');
assert(all(data.MW_g_per_mol>0),'All molecular weights must be positive.');
assert(all(isfinite(data.logP)) && all(isfinite(data.log10_P0_cm_per_s)), ...
    'The regression input contains missing or nonfinite values.');

data.Descriptor_logP_minus_half_log10_MW = ...
    data.logP-0.5*log10(data.MW_g_per_mol);
model=fitlm(data,'log10_P0_cm_per_s ~ Descriptor_logP_minus_half_log10_MW');
data.Predicted_log10_P0_cm_per_s=predict(model,data);
data.Residual_log10_P0=data.log10_P0_cm_per_s-data.Predicted_log10_P0_cm_per_s;
data.Predicted_P0_cm_per_s=10.^data.Predicted_log10_P0_cm_per_s;

intercept=model.Coefficients.Estimate(1);
slope=model.Coefficients.Estimate(2);
ci=coefCI(model,0.05);
rmse=sqrt(mean(data.Residual_log10_P0.^2));
n=height(data);

% Leave-one-out cross-validation using the identical model structure.
looPrediction=nan(n,1);
for i=1:n
    training=data; training(i,:)=[];
    looModel=fitlm(training,'log10_P0_cm_per_s ~ Descriptor_logP_minus_half_log10_MW');
    looPrediction(i)=predict(looModel,data(i,:));
end
data.LOO_Predicted_log10_P0_cm_per_s=looPrediction;
data.LOO_Residual_log10_P0=data.log10_P0_cm_per_s-looPrediction;
q2=corr(data.log10_P0_cm_per_s,looPrediction)^2;
rmseCV=sqrt(mean(data.LOO_Residual_log10_P0.^2));

equation=string(sprintf('P0 = 10^(%.15g %+.15g*(logP - 0.5*log10(MW)))',intercept,slope));
summary=table("Tsinman-based PAMPA-BBB",n,intercept,slope, ...
    ci(1,1),ci(1,2),ci(2,1),ci(2,2),model.Rsquared.Ordinary, ...
    model.Rsquared.Adjusted,rmse,q2,rmseCV,min(data.logP),max(data.logP), ...
    min(data.MW_g_per_mol),max(data.MW_g_per_mol),equation, ...
    'VariableNames',{'Model','N','Intercept','Slope','Intercept_CI95_Lower', ...
    'Intercept_CI95_Upper','Slope_CI95_Lower','Slope_CI95_Upper','R2', ...
    'Adjusted_R2','RMSE_log10','LOO_Q2','LOO_RMSE_log10','Minimum_logP', ...
    'Maximum_logP','Minimum_MW_g_per_mol','Maximum_MW_g_per_mol','Equation'});

summaryFile=fullfile(outputDir,'tsinman_pampa_model_summary.csv');
predictionFile=fullfile(outputDir,'tsinman_pampa_model_predictions.csv');
writetable(summary,summaryFile); writetable(data,predictionFile);

fprintf('\nTsinman-based PAMPA-BBB model (n = %d)\n',n);
fprintf('log10(P0 [cm/s]) = %.6f %+.6f*(logP - 0.5*log10(MW [g/mol]))\n',intercept,slope);
fprintf('R2 = %.6f; RMSE = %.6f log10 units\n',model.Rsquared.Ordinary,rmse);
fprintf('LOO Q2 = %.6f; LOO RMSE = %.6f log10 units\n',q2,rmseCV);
fprintf('Summary saved: %s\nPredictions saved: %s\n',summaryFile,predictionFile);
