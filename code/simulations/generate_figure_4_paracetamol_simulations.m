%% Fit paracetamol BBB permeability and generate the simulations used in Figure 4

clear; clc;
scriptDir = fileparts(mfilename('fullpath'));
repoRoot = fileparts(fileparts(scriptDir));
addpath(fullfile(scriptDir,'models'));
addpath(fullfile(repoRoot,'code','utilities'));

files.tableS1 = fullfile(repoRoot,'data','model_parameters','table_s1_base_parameters.csv');
files.tableS2 = fullfile(repoRoot,'data','model_parameters','table_s2_plasma_pk_parameters.csv');
files.tableS3 = fullfile(repoRoot,'data','model_parameters','table_s3_paracetamol_base_parameters.csv');
files.tableS4 = fullfile(repoRoot,'data','clinical','table_s4_clinical_studies.csv');
files.settings = fullfile(repoRoot,'data','model_parameters','paracetamol_simulation_settings.csv');
files.scenarios = fullfile(repoRoot,'results','tables','figure_3_permeability.csv');
files.clinical = fullfile(repoRoot,'data','clinical','figure_4','Paracetamol_clinical_data.csv');
for fileName = string(struct2cell(files))'
    if ~isfile(fileName), error('Required input was not found: %s',fileName); end
end

requiredFunctions = ["setoptcompphysiol","initdrugdbCNS","drugtemplateCNS", ...
    "obstemplateCNS","plottemplateCNS","FiveBrainobservables","Physiology", ...
    "Individual","loaddrugdata","paracetamol_Wang_plasma_5BRAIN"];
for functionName = requiredFunctions
    if isempty(which(functionName)), error('Required function %s is not on the MATLAB path.',functionName); end
end

outputSimulationDir = fullfile(repoRoot,'results','simulations','main_figures','figure_4');
outputParameterDir = fullfile(repoRoot,'results','model_parameters');
outputTableDir = fullfile(repoRoot,'results','tables');
if ~exist(outputSimulationDir,'dir'), mkdir(outputSimulationDir); end
if ~exist(outputParameterDir,'dir'), mkdir(outputParameterDir); end
if ~exist(outputTableDir,'dir'), mkdir(outputTableDir); end
fitResultFile = fullfile(outputParameterDir,'paracetamol_optimized_permeability.csv');
simulationOutputFile = fullfile(outputSimulationDir,'Paracetamol_simulations.csv');
tableS3OutputFile = fullfile(outputTableDir,'table_s3_paracetamol.csv');

configureCompphysiol();
inputs = loadInputs(files);
indv = buildIndividual(inputs);
fit = optimizePSBBB(indv,inputs.clinicalCSF,inputs.optimization);
fitTable = table(fit.optimizedPS,fit.minimumAAFE,fit.exitflag,fit.iterations, ...
    inputs.optimization.reportedLowerBound,inputs.optimization.reportedUpperBound,false,"spinal CSF", ...
    'VariableNames',{'OptimizedPS_BBB_L_per_h','MinimumAAFE','ExitFlag','Iterations', ...
    'ReportedLowerBound_L_per_h','ReportedUpperBound_L_per_h','BoundsApplied','OptimizationTarget'});
writetable(fitTable,fitResultFile);

tableS3Paracetamol = build_table_s3_paracetamol(files.tableS3,files.tableS1,fitResultFile,files.settings);
writetable(tableS3Paracetamol,tableS3OutputFile);
scenarios = [inputs.scenarios;table("optimized",fit.optimizedPS,'VariableNames',{'Source','PS_BBB_L_per_h'}); ...
    table("perfusion-limited",inputs.perfusionLimitedPS,'VariableNames',{'Source','PS_BBB_L_per_h'})];
simulations = simulateScenarios(indv,scenarios);
writetable(simulations,simulationOutputFile);

fprintf('Optimized PS_BBB: %.15g L/h; minimum AAFE: %.15g\n',fit.optimizedPS,fit.minimumAAFE);
fprintf('Table S3 saved: %s\nSimulations saved: %s\n',tableS3OutputFile,simulationOutputFile);

%% Local workflow functions
function configureCompphysiol()
setoptcompphysiol('DrugDB',@initdrugdbCNS);
setoptcompphysiol('DrugTemplate',@drugtemplateCNS);
setoptcompphysiol('ObservableTemplate',@obstemplateCNS);
setoptcompphysiol('PlotTemplate',@plottemplateCNS);
setoptcompphysiol('DisplayUnits',{'mg/L','L/h'});
end

function inputs = loadInputs(files)
s1 = readtable(files.tableS1,'TextType','string');
s2 = readStrings(files.tableS2); s3 = readStrings(files.tableS3);
s4 = readtable(files.tableS4,'TextType','string'); settings = readStrings(files.settings);
s1v = @(name) s1.Value(s1.Parameter == name);
s2v = @(name) str2double(s2.Paracetamol_value(s2.Parameter == name));
s3v = @(name) str2double(s3.Value(s3.Parameter == name));
sv = @(name) str2double(settings.Value(settings.Parameter == name));
names = ["V_brt","V_brb","endothelialFraction","csfFraction","cranialCSFFraction", ...
    "brainECFFraction","Qcbf","Qproduction","bulkFlowFraction","spinalSinkFraction", ...
    "spinalOutflowFraction","BCSFBtoBBBRatio","BCMtoBBBRatio"];
keys = ["brain_tissue_volume","brain_blood_volume","endothelial_fraction_of_total_brain", ...
    "csf_fraction_of_total_brain","cranial_fraction_of_total_csf","brain_ecf_fraction", ...
    "cerebral_blood_flow","csf_production_rate","bulk_flow_fraction","spinal_sink_fraction", ...
    "spinal_outflow_fraction","bcsfb_to_bbb_surface_area_ratio","bcm_to_bbb_surface_area_ratio"];
for i=1:numel(names), inputs.physiology.(names(i)) = s1v(keys(i)); end
inputs.pk = struct('ka',s2v("ka"),'F',s2v("F"),'Tlag',s2v("Tlag"), ...
    'V1',s2v("Vcntr"),'V2',s2v("Vpone"),'V3',s2v("Vptwo"), ...
    'Q2',s2v("Qcpone"),'Q3',s2v("Qcptwo"),'CL',s2v("CL"));
inputs.drug = struct('MW',s3v("MW"),'BP',s3v("BP"),'fuPlasma',s3v("fu_plasma"), ...
    'fuBrainPlasma',s3v("fu_brain_plasma"),'fuBrainECF',s3v("fu_brain_ecf"), ...
    'fuBrainICF',s3v("fu_brain_icf"),'fuCSF',s3v("fu_csf"), ...
    'fnBrainPlasma',s3v("fn_brain_plasma"),'fnBrainECF',s3v("fn_brain_ecf"), ...
    'fnBrainICF',s3v("fn_brain_icf"),'fnCSF',s3v("fn_csf"));
ivRows = s4(s4.Drug == "Paracetamol" & contains(s4.DosingRegimen,"IV"),:);
assert(numel(unique(ivRows.ModelDoseValue))==1 && numel(unique(ivRows.ModelAdministrationDurationMin))==1, ...
    'Table S4 contains inconsistent paracetamol IV modeling regimens.');
inputs.dose.valueMg = ivRows.ModelDoseValue(1)*1000; inputs.dose.durationMin = ivRows.ModelAdministrationDurationMin(1);
inputs.PSe = sv("PS_ependyma"); inputs.perfusionLimitedPS = sv("perfusion_limited_PS_BBB");
inputs.samplingEnd = sv("sampling_end"); inputs.samplingInterval = sv("sampling_interval");
inputs.optimization.initial = sv("initial_PS_BBB");
inputs.optimization.reportedLowerBound = sv("reported_lower_PS_BBB");
inputs.optimization.reportedUpperBound = sv("reported_upper_PS_BBB");
scenarios = readtable(files.scenarios,'TextType','string');
inputs.scenarios = scenarios(scenarios.Compound == "Paracetamol",{'Source','PS_BBB_L_per_h'});
inputs.scenarios.Source(inputs.scenarios.Source == "in vitro MDCK") = "in vitro cell-based";
clinical = readtable(files.clinical,'TextType','string'); inputs.clinicalCSF = clinical(clinical.site == "scsf",:);
end

function indv = buildIndividual(inputs)
humanphys = Physiology('human35m'); BW = getvalue(humanphys,'BW'); ratio = double(BW/(70*u.kg));
par = parameters('lambda_po',0/u.h,'Tlag',0*u.h,'F',0, ...
    'V1',inputs.pk.V1*ratio*u.L,'V2',inputs.pk.V2*ratio*u.L,'V3',inputs.pk.V3*ratio*u.L, ...
    'Q2',inputs.pk.Q2*ratio*u.L/u.h,'Q3',inputs.pk.Q3*ratio*u.L/u.h,'CL',inputs.pk.CL*ratio^0.75*u.L/u.h, ...
    'PS_BBB',inputs.optimization.initial*u.L/u.h, ...
    'V_brt',inputs.physiology.V_brt*u.L,'V_brb',inputs.physiology.V_brb*u.L, ...
    'endothelial_fraction',inputs.physiology.endothelialFraction,'csf_fraction',inputs.physiology.csfFraction, ...
    'cranial_csf_fraction',inputs.physiology.cranialCSFFraction,'brain_ecf_fraction',inputs.physiology.brainECFFraction, ...
    'Q_cbf',inputs.physiology.Qcbf*u.L/u.h,'Q_production',inputs.physiology.Qproduction*u.L/u.h, ...
    'bulk_flow_fraction',inputs.physiology.bulkFlowFraction,'spinal_sink_fraction',inputs.physiology.spinalSinkFraction, ...
    'spinal_outflow_fraction',inputs.physiology.spinalOutflowFraction, ...
    'BCSFB_to_BBB_ratio',inputs.physiology.BCSFBtoBBBRatio,'BCM_to_BBB_ratio',inputs.physiology.BCMtoBBBRatio, ...
    'PS_e',inputs.PSe*u.L/u.h,'MW',inputs.drug.MW*u.g/u.mol,'BP',inputs.drug.BP, ...
    'fu_plasma',inputs.drug.fuPlasma,'fu_brain_plasma',inputs.drug.fuBrainPlasma, ...
    'fu_brain_ecf',inputs.drug.fuBrainECF,'fu_brain_icf',inputs.drug.fuBrainICF,'fu_csf',inputs.drug.fuCSF, ...
    'fn_brain_plasma',inputs.drug.fnBrainPlasma,'fn_brain_ecf',inputs.drug.fnBrainECF, ...
    'fn_brain_icf',inputs.drug.fnBrainICF,'fn_csf',inputs.drug.fnCSF);
indv = Individual('Virtual',1); indv.name='Human'; indv.physiology=humanphys;
indv.dosing=Infusion('Paracetamol',0*u.h,inputs.dose.valueMg*u.mg,inputs.dose.durationMin*u.min,'iv');
indv.drugdata=loaddrugdata('Paracetamol','species','human');
indv.sampling=Sampling((0:inputs.samplingInterval:inputs.samplingEnd)*u.h,FiveBrainobservables);
indv.model=paracetamol_Wang_plasma_5BRAIN; indv.model.par=par;
end

function fit = optimizePSBBB(indv,clinical,settings)
objective = @(logPS) aafe(logPS,indv,clinical.time,clinical.conc,1e-6);
options = optimset('Display','iter','TolFun',1e-6,'TolX',1e-6);
[logPS,fval,exitflag,output] = fminsearch(objective,log(settings.initial),options);
fit = struct('optimizedPS',exp(logPS),'minimumAAFE',fval,'exitflag',exitflag,'iterations',output.iterations);
end

function value = aafe(logPS,indv,timeData,concentrationData,epsilon)
local = indv;
try
    local.model.par.PS_BBB=exp(logPS)*u.L/u.h; initialize(local); simulate(local);
    t=double(local.observation.Time(8:8:end))/3600; c=double(local.observation.Value(8:8:end))*1000;
    pred=interp1(t,c,timeData,'linear','extrap'); value=exp(mean(abs(log(pred+epsilon)-log(concentrationData+epsilon))));
catch
    value=1e6;
end
end

function simulations = simulateScenarios(indv,scenarios)
initialize(indv); simulate(indv); tAll=double(indv.observation.Time)/3600; v=double(indv.observation.Value);
t=tAll(1:8:end); n=numel(t);
plasma=table(repmat("plasma",n,1),NaN(n,1),repmat("",n,1),t,v(1:8:end)*1000, ...
    'VariableNames',{'site','PS','source','time','conc'});
sites=["brain ECF","scsf"]; offsets=[4,8]; out=cell(height(scenarios)*2,1); k=0;
for j=1:2
    for i=1:height(scenarios)
        indv.model.par.PS_BBB=scenarios.PS_BBB_L_per_h(i)*u.L/u.h; initialize(indv); simulate(indv); v=double(indv.observation.Value); k=k+1;
        out{k}=table(repmat(sites(j),n,1),repmat(scenarios.PS_BBB_L_per_h(i),n,1), ...
            repmat(scenarios.Source(i),n,1),t,v(offsets(j):8:end)*1000, ...
            'VariableNames',{'site','PS','source','time','conc'});
    end
end
simulations=[plasma;vertcat(out{:})];
end

function data = readStrings(file)
opts=detectImportOptions(file); opts=setvartype(opts,opts.VariableNames,'string'); data=readtable(file,opts);
end
