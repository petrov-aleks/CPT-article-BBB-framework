%% Fit indomethacin BBB permeability and generate the simulations used in Figure 4

clear; clc;

scriptDir = fileparts(mfilename('fullpath'));
repoRoot = fileparts(fileparts(scriptDir));
addpath(fullfile(scriptDir,'models'));
addpath(fullfile(repoRoot,'code','utilities'));

files.tableS1 = fullfile(repoRoot,'data','model_parameters','table_s1_base_parameters.csv');
files.tableS2 = fullfile(repoRoot,'data','model_parameters','table_s2_plasma_pk_parameters.csv');
files.tableS3 = fullfile(repoRoot,'data','model_parameters','table_s3_indomethacin_base_parameters.csv');
files.tableS4 = fullfile(repoRoot,'data','clinical','table_s4_clinical_studies.csv');
files.settings = fullfile(repoRoot,'data','model_parameters','indomethacin_simulation_settings.csv');
files.scenarios = fullfile(repoRoot,'results','tables','figure_3_permeability.csv');
files.clinical = fullfile(repoRoot,'data','clinical','figure_4','Indomethacin_clinical_data.csv');

fileNames = struct2cell(files);
for i = 1:numel(fileNames)
    if ~isfile(fileNames{i}), error('Required input was not found: %s',fileNames{i}); end
end

requiredFunctions = ["setoptcompphysiol","initdrugdbCNS","drugtemplateCNS", ...
    "obstemplateCNS","plottemplateCNS","FiveBrainobservables","Physiology", ...
    "Individual","loaddrugdata","indometacin_Saleh_plasma_5BRAIN","lsqnonlin"];
for functionName = requiredFunctions
    if isempty(which(functionName)), error('Required function %s is not on the MATLAB path.',functionName); end
end

outputSimulationDir = fullfile(repoRoot,'results','simulations','main_figures','figure_4');
outputParameterDir = fullfile(repoRoot,'results','model_parameters');
outputTableDir = fullfile(repoRoot,'results','tables');
if ~exist(outputSimulationDir,'dir'), mkdir(outputSimulationDir); end
if ~exist(outputParameterDir,'dir'), mkdir(outputParameterDir); end
if ~exist(outputTableDir,'dir'), mkdir(outputTableDir); end

fitResultFile = fullfile(outputParameterDir,'indomethacin_optimized_permeability.csv');
simulationOutputFile = fullfile(outputSimulationDir,'Indomethacin_simulations.csv');
tableS3OutputFile = fullfile(outputTableDir,'table_s3_indomethacin.csv');

configureCompphysiol();
inputs = loadIndomethacinInputs(files);
indv = buildIndomethacinIndividual(inputs);

fit = optimizePSBBB(indv,inputs.clinicalCSF,inputs.optimization);
fitTable = table(fit.optimizedPS,fit.resnorm,fit.exitflag, ...
    inputs.optimization.lowerBound,inputs.optimization.upperBound,"spinal CSF", ...
    'VariableNames',{'OptimizedPS_BBB_L_per_h','ResidualNorm','ExitFlag', ...
    'LowerBound_L_per_h','UpperBound_L_per_h','OptimizationTarget'});
writetable(fitTable,fitResultFile);

tableS3Indomethacin = build_table_s3_indomethacin(files.tableS3,files.tableS1, ...
    fitResultFile,files.settings);
writetable(tableS3Indomethacin,tableS3OutputFile);

scenarios = [inputs.scenarios; ...
    table("optimized",fit.optimizedPS,'VariableNames',{'Source','PS_BBB_L_per_h'}); ...
    table("perfusion-limited",inputs.perfusionLimitedPS,'VariableNames',{'Source','PS_BBB_L_per_h'})];
indomethacinSimulations = simulateScenarios(indv,scenarios);
writetable(indomethacinSimulations,simulationOutputFile);


fprintf('Optimized PS_BBB: %.15g L/h\n',fit.optimizedPS);
fprintf('Table S3 saved: %s\n',tableS3OutputFile);
fprintf('Simulations saved: %s\n',simulationOutputFile);

%% Local workflow functions

function configureCompphysiol()
setoptcompphysiol('DrugDB',@initdrugdbCNS);
setoptcompphysiol('DrugTemplate',@drugtemplateCNS);
setoptcompphysiol('ObservableTemplate',@obstemplateCNS);
setoptcompphysiol('PlotTemplate',@plottemplateCNS);
setoptcompphysiol('DisplayUnits',{'ng/mL','L/h'});
end

function inputs = loadIndomethacinInputs(files)
s1 = readtable(files.tableS1,'TextType','string');
s2 = readAllAsString(files.tableS2);
s3 = readAllAsString(files.tableS3);
s4 = readtable(files.tableS4,'TextType','string');
settings = readAllAsString(files.settings);
s1value = @(name) s1.Value(s1.Parameter == name);
s2value = @(name) str2double(s2.Indomethacin_value(s2.Parameter == name));
s3value = @(name) str2double(s3.Value(s3.Parameter == name));
setting = @(name) str2double(settings.Value(settings.Parameter == name));

inputs.physiology.V_brt = s1value("brain_tissue_volume")*u.L;
inputs.physiology.V_brb = s1value("brain_blood_volume")*u.L;
inputs.physiology.endothelialFraction = s1value("endothelial_fraction_of_total_brain");
inputs.physiology.csfFraction = s1value("csf_fraction_of_total_brain");
inputs.physiology.cranialCSFFraction = s1value("cranial_fraction_of_total_csf");
inputs.physiology.brainECFFraction = s1value("brain_ecf_fraction");
inputs.physiology.Qcbf = s1value("cerebral_blood_flow")*u.L/u.h;
inputs.physiology.Qproduction = s1value("csf_production_rate")*u.L/u.h;
inputs.physiology.bulkFlowFraction = s1value("bulk_flow_fraction");
inputs.physiology.spinalSinkFraction = s1value("spinal_sink_fraction");
inputs.physiology.spinalOutflowFraction = s1value("spinal_outflow_fraction");
inputs.physiology.BCSFBtoBBBRatio = s1value("bcsfb_to_bbb_surface_area_ratio");
inputs.physiology.BCMtoBBBRatio = s1value("bcm_to_bbb_surface_area_ratio");

inputs.pk.ka = s2value("ka");
inputs.pk.F = s2value("F");
inputs.pk.V1 = s2value("Vcntr");
inputs.pk.V2 = s2value("Vpone");
inputs.pk.CL1 = s2value("CL");
inputs.pk.CL2 = s2value("Qcpone");
inputs.pk.Tlag = s2value("Tlag");

inputs.drug.MW = s3value("MW")*u.g/u.mol;
inputs.drug.BP = s3value("BP");
inputs.drug.fuPlasma = s3value("fu_plasma");
inputs.drug.fuBrainPlasma = s3value("fu_brain_plasma");
inputs.drug.fuBrainECF = s3value("fu_brain_ecf");
inputs.drug.fuBrainICF = s3value("fu_brain_icf");
qAlbumin = s3value("albumin_csf_to_plasma_ratio");
inputs.drug.fuCSF = 1/(1+qAlbumin*(1/inputs.drug.fuPlasma-1));
fractionNeutral = @(pH) 1/(1+10^(pH-s3value("pKa")));
inputs.drug.fnBrainPlasma = fractionNeutral(s3value("pH_brain_plasma"));
inputs.drug.fnBrainECF = fractionNeutral(s3value("pH_brain_ecf"));
inputs.drug.fnBrainICF = fractionNeutral(s3value("pH_brain_icf"));
inputs.drug.fnCSF = fractionNeutral(s3value("pH_csf"));

indomethacinRows = s4(s4.Drug == "Indomethacin",:);
assert(numel(unique(indomethacinRows.ModelDoseValue)) == 1 && ...
    numel(unique(indomethacinRows.ModelDoseUnit)) == 1 && ...
    numel(unique(indomethacinRows.ModelAdministrationDurationMin)) == 1, ...
    'Table S4 contains inconsistent indomethacin modeling regimens.');
inputs.dose.valueMg = indomethacinRows.ModelDoseValue(1);

inputs.PSe = setting("PS_ependyma")*u.L/u.h;
inputs.perfusionLimitedPS = setting("perfusion_limited_PS_BBB");
inputs.samplingEnd = setting("sampling_end");
inputs.samplingInterval = setting("sampling_interval");
inputs.optimization.initial = setting("initial_PS_BBB");
inputs.optimization.lowerBound = setting("lower_PS_BBB");
inputs.optimization.upperBound = setting("upper_PS_BBB");

scenarios = readtable(files.scenarios,'TextType','string');
scenarios = scenarios(scenarios.Compound == "Indomethacin",{'Source','PS_BBB_L_per_h'});
scenarios.Source(scenarios.Source == "in vitro Caco-2") = "in vitro cell-based";
inputs.scenarios = scenarios;

clinical = readtable(files.clinical,'TextType','string');
inputs.clinicalCSF = clinical(clinical.site == "scsf",:);
end

function indv = buildIndomethacinIndividual(inputs)
humanphys = Physiology('human35m');
par = parameters( ...
    'V1',inputs.pk.V1*u.L, ...
    'V2',inputs.pk.V2*u.L, ...
    'CL1',inputs.pk.CL1*u.L/u.h,'CL2',inputs.pk.CL2*u.L/u.h, ...
    'Tlag',inputs.pk.Tlag*u.h, ...
    'lambda_po',inputs.pk.ka/u.h,'F',inputs.pk.F, ...
    'PS_BBB',inputs.optimization.initial*u.L/u.h, ...
    'V_brt',inputs.physiology.V_brt,'V_brb',inputs.physiology.V_brb, ...
    'endothelial_fraction',inputs.physiology.endothelialFraction, ...
    'csf_fraction',inputs.physiology.csfFraction, ...
    'cranial_csf_fraction',inputs.physiology.cranialCSFFraction, ...
    'brain_ecf_fraction',inputs.physiology.brainECFFraction, ...
    'Q_cbf',inputs.physiology.Qcbf,'Q_production',inputs.physiology.Qproduction, ...
    'bulk_flow_fraction',inputs.physiology.bulkFlowFraction, ...
    'spinal_sink_fraction',inputs.physiology.spinalSinkFraction, ...
    'spinal_outflow_fraction',inputs.physiology.spinalOutflowFraction, ...
    'BCSFB_to_BBB_ratio',inputs.physiology.BCSFBtoBBBRatio, ...
    'BCM_to_BBB_ratio',inputs.physiology.BCMtoBBBRatio,'PS_e',inputs.PSe, ...
    'MW',inputs.drug.MW,'BP',inputs.drug.BP,'fu_plasma',inputs.drug.fuPlasma, ...
    'fu_brain_plasma',inputs.drug.fuBrainPlasma,'fu_brain_ecf',inputs.drug.fuBrainECF, ...
    'fu_brain_icf',inputs.drug.fuBrainICF,'fu_csf',inputs.drug.fuCSF, ...
    'fn_brain_plasma',inputs.drug.fnBrainPlasma,'fn_brain_ecf',inputs.drug.fnBrainECF, ...
    'fn_brain_icf',inputs.drug.fnBrainICF,'fn_csf',inputs.drug.fnCSF);

indv = Individual('Virtual',1);
indv.name = 'Human';
indv.physiology = humanphys;
% Oral is the toolbox implementation of first-order extravascular input;
% the clinical administration route represented here is intramuscular.
indv.dosing = Oral('Indometacin',0*u.h,inputs.dose.valueMg*u.mg);
indv.drugdata = loaddrugdata('Indometacin','species','human');
indv.sampling = Sampling((0:inputs.samplingInterval:inputs.samplingEnd)*u.h,FiveBrainobservables);
indv.model = indometacin_Saleh_plasma_5BRAIN;
indv.model.par = par;
end

function fit = optimizePSBBB(indv,clinicalCSF,settings)
objective = @(logPS) spinalCSFLogResiduals(logPS,indv,clinicalCSF.time,clinicalCSF.conc,1e-6);
options = optimoptions('lsqnonlin','Display','iter','MaxFunctionEvaluations',200, ...
    'FunctionTolerance',1e-6);
[logPS,resnorm,~,exitflag,output] = lsqnonlin(objective,log(settings.initial), ...
    log(settings.lowerBound),log(settings.upperBound),options);
fit.optimizedPS = exp(logPS);
fit.resnorm = resnorm;
fit.exitflag = exitflag;
fit.output = output;
end

function residuals = spinalCSFLogResiduals(logPS,indv,timeData,concentrationData,epsilon)
localIndividual = indv;
try
    localIndividual.model.par.PS_BBB = exp(logPS)*u.L/u.h;
    initialize(localIndividual); simulate(localIndividual);
    simulationTime = double(localIndividual.observation.Time(8:8:end))/3600;
    spinalCSF = double(localIndividual.observation.Value(8:8:end))*1e6;
    predictions = interp1(simulationTime,spinalCSF,timeData,'linear','extrap');
    residuals = log(concentrationData+epsilon)-log(predictions+epsilon);
catch
    residuals = 1e3*ones(size(concentrationData));
end
end

function simulations = simulateScenarios(indv,scenarios)
initialize(indv); simulate(indv);
allTimes = double(indv.observation.Time)/3600;
allValues = double(indv.observation.Value);
simulationTime = allTimes(1:8:end);
nTime = numel(simulationTime);
plasma = table(repmat("plasma",nTime,1),NaN(nTime,1),repmat("",nTime,1), ...
    simulationTime,allValues(1:8:end)*1e6,'VariableNames',{'site','PS','source','time','conc'});
sites = ["brain ECF","scsf"];
offsets = [4,8];
tables = cell(height(scenarios)*2,1);
index = 0;
for siteIndex = 1:2
    for scenarioIndex = 1:height(scenarios)
        indv.model.par.PS_BBB = scenarios.PS_BBB_L_per_h(scenarioIndex)*u.L/u.h;
        initialize(indv); simulate(indv);
        values = double(indv.observation.Value);
        index = index+1;
        tables{index} = table(repmat(sites(siteIndex),nTime,1), ...
            repmat(scenarios.PS_BBB_L_per_h(scenarioIndex),nTime,1), ...
            repmat(scenarios.Source(scenarioIndex),nTime,1),simulationTime, ...
            values(offsets(siteIndex):8:end)*1e6,'VariableNames',{'site','PS','source','time','conc'});
    end
end
simulations = [plasma;vertcat(tables{:})];
end

function data = readAllAsString(file)
options = detectImportOptions(file);
options = setvartype(options,options.VariableNames,'string');
data = readtable(file,options);
end
