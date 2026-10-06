%% Generate Table S5: published and Tsinman-based in silico models

clear; clc;
scriptDir=fileparts(mfilename('fullpath'));
repoRoot=fileparts(fileparts(scriptDir));
modelFile=fullfile(repoRoot,'results','model_development','tsinman_pampa_model_summary.csv');
outputDir=fullfile(repoRoot,'results','tables');
if ~isfile(modelFile)
    error('Run code/analysis/fit_tsinman_pampa_model.m before generating Table S5.');
end
if ~exist(outputDir,'dir'), mkdir(outputDir); end

fit=readtable(modelFile,'TextType','string');
assert(height(fit)==1,'Expected one fitted Tsinman model summary.');

% Published models are entered directly from their source publications.
category=["In silico perfusion models";"In silico perfusion models"; ...
    "In silico PAMPA models";"In silico PAMPA models"];
model=["Smith et al. 2024";"Luco et al. 2006"; ...
    "Tsinman-based model";"Grumetto et al. 2016"];
equation=[ ...
    "PS0 = 10^(-1.24 + 0.887*(logP - 0.5*log10(MW)))"; ...
    "PS0 = 10^(-2.06 + 0.448*logP - 0.366*(MW/100))"; ...
    "P0 = 10^("+compose('%.2f',fit.Intercept)+" + "+compose('%.2f',fit.Slope)+"*(logP - 0.5*log10(MW)))"; ...
    "P0 = 10^(-6.210 + 0.939*logP)"];
r2=[0.95;0.80;fit.R2;0.78];
n=[78;22;fit.N;36];
descriptor=["logP - 0.5*log10(MW)";"logP, MW"; ...
    "logP - 0.5*log10(MW)";"logP"];
minimumLogP=[-4.4;-3.3;fit.Minimum_logP;-0.8];
maximumLogP=[7.57;4.82;fit.Maximum_logP;5.2];
minimumMW=[18;75;fit.Minimum_MW_g_per_mol;76];
maximumMW=[5000;515;fit.Maximum_MW_g_per_mol;543];
reference=["Smith et al. 2024";"Luco et al. 2006"; ...
    "Tsinman et al. 2011";"Grumetto et al. 2016"];
valueOrigin=["Published equation and statistics";"Published equation and statistics"; ...
    "Equation and statistics fitted in the present work";"Published equation and statistics"];

% Machine-readable coefficients used by downstream prediction scripts. The
% exponent is evaluated as
%   Intercept + LogP_Coefficient*logP ...
%             + Log10_MW_Coefficient*log10(MW) ...
%             + MW_over_100_Coefficient*(MW/100).
% Four decimal places are retained for the fitted Tsinman-based model, matching
% the original prediction analysis; Equation is the manuscript-facing form.
tsinmanIntercept=round(fit.Intercept,4);
tsinmanSlope=round(fit.Slope,4);
intercept=[-1.24;-2.06;tsinmanIntercept;-6.210];
logPCoefficient=[0.887;0.448;tsinmanSlope;0.939];
log10MWCoefficient=[-0.5*0.887;0;-0.5*tsinmanSlope;0];
mwOver100Coefficient=[0;-0.366;0;0];

tableS5=table(category,model,equation,r2,n,descriptor,minimumLogP,maximumLogP, ...
    minimumMW,maximumMW,reference,valueOrigin,intercept,logPCoefficient, ...
    log10MWCoefficient,mwOver100Coefficient, ...
    'VariableNames',{'Category','Model','Equation','R2','N','Descriptor', ...
    'Minimum_logP','Maximum_logP','Minimum_MW_g_per_mol','Maximum_MW_g_per_mol', ...
    'Reference','ValueOrigin','Intercept','LogP_Coefficient', ...
    'Log10_MW_Coefficient','MW_over_100_Coefficient'});
csvFile=fullfile(outputDir,'table_s5.csv');
writetable(tableS5,csvFile);

disp(tableS5)
fprintf('Saved %s\n',csvFile);
