%% Fit ethanol BBB permeability and generate the simulations used in Figure 4

clear; clc;

scriptDir = fileparts(mfilename('fullpath'));
repoRoot = fileparts(fileparts(scriptDir));
addpath(fullfile(scriptDir,'models'));
addpath(fullfile(repoRoot,'code','utilities'));

files.tableS1 = fullfile(repoRoot,'data','model_parameters','table_s1_base_parameters.csv');
files.tableS2 = fullfile(repoRoot,'data','model_parameters','table_s2_plasma_pk_parameters.csv');
files.tableS3 = fullfile(repoRoot,'data','model_parameters','table_s3_ethanol_base_parameters.csv');
files.tableS4 = fullfile(repoRoot,'data','clinical','table_s4_clinical_studies.csv');
files.settings = fullfile(repoRoot,'data','model_parameters','ethanol_simulation_settings.csv');
files.scenarios = fullfile(repoRoot,'results','tables','figure_3_permeability.csv');
files.clinical = fullfile(repoRoot,'data','clinical','figure_4','Ethanol_clinical_data.csv');

fileNames = struct2cell(files);
for i = 1:numel(fileNames)
    if ~isfile(fileNames{i}), error('Required input was not found: %s',fileNames{i}); end
end

requiredFunctions = ["setoptcompphysiol","initdrugdbCNS","drugtemplateCNS", ...
    "obstemplateCNS","plottemplateCNS","FiveBrainobservables","Physiology", ...
    "Individual","loaddrugdata","ethanol_plasma_5BRAIN","lsqnonlin"];
for functionName = requiredFunctions
    if isempty(which(functionName)), error('Required function %s is not on the MATLAB path.',functionName); end
end

outputSimulationDir = fullfile(repoRoot,'results','simulations','main_figures','figure_4');
outputParameterDir = fullfile(repoRoot,'results','model_parameters');
outputTableDir = fullfile(repoRoot,'results','tables');
if ~exist(outputSimulationDir,'dir'), mkdir(outputSimulationDir); end
if ~exist(outputParameterDir,'dir'), mkdir(outputParameterDir); end
if ~exist(outputTableDir,'dir'), mkdir(outputTableDir); end

fitResultFile = fullfile(outputParameterDir,'ethanol_optimized_permeability.csv');
simulationOutputFile = fullfile(outputSimulationDir,'Ethanol_simulations.csv');
tableS3OutputFile = fullfile(outputTableDir,'table_s3_ethanol.csv');

configureCompphysiol();
inputs = loadEthanolInputs(files);
indv = buildEthanolIndividual(inputs);

fit = optimizePSBBB(indv,inputs.clinicalBrain,inputs.optimization);
fitTable = table(fit.optimizedPS,fit.resnorm,fit.exitflag, ...
    inputs.optimization.lowerBound,inputs.optimization.upperBound,"brain ECF", ...
    'VariableNames',{'OptimizedPS_BBB_L_per_h','ResidualNorm','ExitFlag', ...
    'LowerBound_L_per_h','UpperBound_L_per_h','OptimizationTarget'});
writetable(fitTable,fitResultFile);

tableS3Ethanol = build_table_s3_ethanol(files.tableS3,files.tableS1, ...
    fitResultFile,files.settings);
writetable(tableS3Ethanol,tableS3OutputFile);

scenarios = [inputs.scenarios; ...
    table("optimized",fit.optimizedPS,'VariableNames',{'Source','PS_BBB_L_per_h'}); ...
    table("perfusion-limited",inputs.perfusionLimitedPS,'VariableNames',{'Source','PS_BBB_L_per_h'})];
ethanolSimulations = simulateScenarios(indv,scenarios);
writetable(ethanolSimulations,simulationOutputFile);


fprintf('Optimized PS_BBB: %.15g L/h\n',fit.optimizedPS);
fprintf('Table S3 saved: %s\n',tableS3OutputFile);
fprintf('Simulations saved: %s\n',simulationOutputFile);

%% Local workflow functions

function configureCompphysiol()
setoptcompphysiol('DrugDB',@initdrugdbCNS);
setoptcompphysiol('DrugTemplate',@drugtemplateCNS);
setoptcompphysiol('ObservableTemplate',@obstemplateCNS);
setoptcompphysiol('PlotTemplate',@plottemplateCNS);
setoptcompphysiol('DisplayUnits',{'g/L','L/h'});
end

function inputs = loadEthanolInputs(files)
s1 = readtable(files.tableS1,'TextType','string');
s2 = readAllAsString(files.tableS2);
s3 = readAllAsString(files.tableS3);
s4 = readtable(files.tableS4,'TextType','string');
settings = readAllAsString(files.settings);
s1value = @(name) s1.Value(s1.Parameter == name);
s2value = @(name) str2double(s2.Ethanol_value(s2.Parameter == name));
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
inputs.pk.Q = s2value("Qcpone");
inputs.pk.Vmax = s2value("Vmax");
inputs.pk.Km = s2value("KM");
inputs.pk.foodEffectVmax = s2value("FEVmax");
inputs.pk.BWExponentV2 = s2value("BWVpone");
inputs.pk.BWReference = s2value("BWref");

inputs.drug.MW = s3value("MW")*u.g/u.mol;
inputs.drug.BP = 1/s3value("plasma_to_blood_ratio");
inputs.drug.fuPlasma = s3value("fu_plasma");
inputs.drug.fuBrainPlasma = s3value("fu_brain_plasma");
inputs.drug.fuBrainECF = s3value("fu_brain_ecf");
inputs.drug.fuBrainICF = s3value("fu_brain_icf");
inputs.drug.fuCSF = s3value("fu_csf");
inputs.drug.fnBrainPlasma = s3value("fn_brain_plasma");
inputs.drug.fnBrainECF = s3value("fn_brain_ecf");
inputs.drug.fnBrainICF = s3value("fn_brain_icf");
inputs.drug.fnCSF = s3value("fn_csf");

ethanolRows = s4(s4.Drug == "Ethanol",:);
assert(numel(unique(ethanolRows.ModelDoseValue)) == 1 && ...
    numel(unique(ethanolRows.ModelDoseUnit)) == 1 && ...
    numel(unique(ethanolRows.ModelAdministrationDurationMin)) == 1, ...
    'Table S4 contains inconsistent ethanol modeling regimens.');
inputs.dose.valuePerKg = ethanolRows.ModelDoseValue(1);
inputs.dose.durationMin = ethanolRows.ModelAdministrationDurationMin(1);

inputs.PSe = setting("PS_ependyma")*u.L/u.h;
inputs.perfusionLimitedPS = setting("perfusion_limited_PS_BBB");
inputs.samplingEnd = setting("sampling_end");
inputs.samplingInterval = setting("sampling_interval");
inputs.optimization.initial = setting("initial_PS_BBB");
inputs.optimization.lowerBound = setting("minimum_PS_BBB");

scenarios = readtable(files.scenarios,'TextType','string');
scenarios = scenarios(scenarios.Compound == "Ethanol",{'Source','PS_BBB_L_per_h'});
scenarios.Source(scenarios.Source == "ex vivo UC") = "in vitro cell-based";
inputs.scenarios = scenarios;
inputs.optimization.upperBound = max(scenarios.PS_BBB_L_per_h);

clinical = readtable(files.clinical,'TextType','string');
inputs.clinicalBrain = clinical(clinical.site == "brain mass",:);
end

function indv = buildEthanolIndividual(inputs)
humanphys = Physiology('human35m');
BW = getvalue(humanphys,'BW');
par = parameters( ...
    'V1',inputs.pk.V1*u.L, ...
    'V2',inputs.pk.V2*(double(BW)/inputs.pk.BWReference)^inputs.pk.BWExponentV2*u.L, ...
    'Vmax',inputs.pk.Vmax*inputs.pk.foodEffectVmax*u.g/u.h, ...
    'Km',inputs.pk.Km*u.g/u.L,'Q',inputs.pk.Q*u.L/u.h, ...
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

dose = inputs.dose.valuePerKg*u.g/u.kg*BW;
indv = Individual('Virtual',1);
indv.name = 'Human';
indv.physiology = humanphys;
indv.dosing = Infusion('Ethanol',0*u.h,dose,inputs.dose.durationMin*u.min,'drink');
indv.drugdata = loaddrugdata('Ethanol','species','human');
indv.sampling = Sampling((0:inputs.samplingInterval:inputs.samplingEnd)*u.h,FiveBrainobservables);
indv.model = ethanol_plasma_5BRAIN;
indv.model.par = par;
end

function fit = optimizePSBBB(indv,clinicalBrain,settings)
objective = @(logPS) brainECFLogResiduals(logPS,indv,clinicalBrain.time,clinicalBrain.conc,1e-6);
options = optimoptions('lsqnonlin','Display','iter','MaxFunctionEvaluations',200,'FunctionTolerance',1e-3);
[logPS,resnorm,~,exitflag,output] = lsqnonlin(objective,log(settings.initial), ...
    log(settings.lowerBound),log(settings.upperBound),options);
fit.optimizedPS = exp(logPS);
fit.resnorm = resnorm;
fit.exitflag = exitflag;
fit.output = output;
end

function residuals = brainECFLogResiduals(logPS,indv,timeData,concentrationData,epsilon)
localIndividual = indv;
try
    localIndividual.model.par.PS_BBB = exp(logPS)*u.L/u.h;
    initialize(localIndividual); simulate(localIndividual);
    simulationTime = double(localIndividual.observation.Time(4:8:end))/3600;
    brainECF = double(localIndividual.observation.Value(4:8:end));
    predictions = interp1(simulationTime,brainECF,timeData,'linear','extrap');
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
    simulationTime,allValues(1:8:end),'VariableNames',{'site','PS','source','time','conc'});
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
            values(offsets(siteIndex):8:end),'VariableNames',{'site','PS','source','time','conc'});
    end
end
simulations = [plasma;vertcat(tables{:})];
end

function data = readAllAsString(file)
options = detectImportOptions(file);
options = setvartype(options,options.VariableNames,'string');
data = readtable(file,options);
end
