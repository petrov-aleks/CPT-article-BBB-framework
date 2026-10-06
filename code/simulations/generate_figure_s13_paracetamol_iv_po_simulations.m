%% Generate paracetamol IV and PO simulations used in Figure S13

clear; clc;
scriptDir=fileparts(mfilename('fullpath')); repoRoot=fileparts(fileparts(scriptDir));
addpath(fullfile(scriptDir,'models'));
setoptcompphysiol('DrugDB',@initdrugdbCNS); setoptcompphysiol('DrugTemplate',@drugtemplateCNS);
setoptcompphysiol('ObservableTemplate',@obstemplateCNS); setoptcompphysiol('PlotTemplate',@plottemplateCNS);
setoptcompphysiol('DisplayUnits',{'g/L','L/h'});

files.s1=fullfile(repoRoot,'data','model_parameters','table_s1_base_parameters.csv');
files.s2=fullfile(repoRoot,'data','model_parameters','table_s2_plasma_pk_parameters.csv');
files.s3=fullfile(repoRoot,'results','tables','table_s3_paracetamol.csv');
files.s4=fullfile(repoRoot,'data','clinical','table_s4_clinical_studies.csv');
files.settings=fullfile(repoRoot,'data','model_parameters','paracetamol_simulation_settings.csv');
files.scenarios=fullfile(repoRoot,'results','tables','figure_3_permeability.csv');
for f=string(struct2cell(files))', if ~isfile(f), error('Required input was not found: %s',f); end, end

inputs=loadInputs(files); outputDir=fullfile(repoRoot,'results','simulations', ...
    'supplementary_figures','figure_s13'); if ~exist(outputDir,'dir'), mkdir(outputDir); end
for route=["IV","PO"]
    individual=buildIndividual(inputs,route);
    simulations=simulateScenarios(individual,inputs.scenarios);
    outputFile=fullfile(outputDir,"Paracetamol_"+route+"_simulations.csv");
    writetable(simulations,outputFile); fprintf('Saved %s\n',outputFile);
end

function in=loadInputs(files)
s1=readtable(files.s1,'TextType','string'); s2=readStrings(files.s2); s3=readStrings(files.s3);
s4=readtable(files.s4,'TextType','string'); settings=readStrings(files.settings);
v1=@(p)s1.Value(s1.Parameter==p); v2=@(p)str2double(s2.Paracetamol_value(s2.Parameter==p));
v3=@(p)exactValue(s3,p); vs=@(p)str2double(settings.Value(settings.Parameter==p));
in.phys=struct('Vbrt',v1("brain_tissue_volume"),'Vbrb',v1("brain_blood_volume"), ...
 'endo',v1("endothelial_fraction_of_total_brain"),'csf',v1("csf_fraction_of_total_brain"), ...
 'cranial',v1("cranial_fraction_of_total_csf"),'ecf',v1("brain_ecf_fraction"), ...
 'Qcbf',v1("cerebral_blood_flow"),'Qprod',v1("csf_production_rate"), ...
 'bulk',v1("bulk_flow_fraction"),'sink',v1("spinal_sink_fraction"), ...
 'out',v1("spinal_outflow_fraction"),'bcsfb',v1("bcsfb_to_bbb_surface_area_ratio"), ...
 'bcm',v1("bcm_to_bbb_surface_area_ratio"));
in.pk=struct('ka',v2("ka"),'F',v2("F"),'Tlag',v2("Tlag"),'V1',v2("Vcntr"), ...
 'V2',v2("Vpone"),'V3',v2("Vptwo"),'Q2',v2("Qcpone"),'Q3',v2("Qcptwo"), ...
 'CL',v2("CL"),'BWref',v2("BWref"),'eV1',v2("BWVcntr"),'eV2',v2("BWVpone"), ...
 'eV3',v2("BWVptwo"),'eQ2',v2("BWQcpone"),'eQ3',v2("BWQcptwo"),'eCL',v2("BWCL"));
in.drug=struct('MW',v3("MW"),'BP',v3("BP"),'fuP',v3("fu_plasma"), ...
 'fuBP',v3("fu_brain_plasma"),'fuE',v3("fu_brain_ECF"),'fuI',v3("fu_brain_ICF"), ...
 'fuC',v3("fu_CSF"),'fnBP',v3("fn_brain_plasma"),'fnE',v3("fn_brain_ECF"), ...
 'fnI',v3("fn_brain_ICF"),'fnC',v3("fn_CSF"),'PS',v3("PS_BBB"),'PSE',v3("PS_E"));
iv=s4(s4.Drug=="Paracetamol" & contains(s4.DosingRegimen,"IV"),:);
po=s4(s4.Drug=="Paracetamol" & contains(s4.DosingRegimen,"PO"),:);
assert(isscalar(unique(iv.ModelDoseValue)) && isscalar(unique(iv.ModelAdministrationDurationMin)));
assert(height(po)==1 && po.ModelDoseUnit=="g");
in.doseIV=iv.ModelDoseValue(1); in.durationIV=iv.ModelAdministrationDurationMin(1);
in.dosePO=po.ModelDoseValue(1); in.dt=vs("sampling_interval"); in.tend=vs("sampling_end");
sc=readtable(files.scenarios,'TextType','string'); sc=sc(sc.Compound=="Paracetamol",{'Source','PS_BBB_L_per_h'});
sc.Source(sc.Source=="in vitro MDCK")="in vitro cell-based";
in.scenarios=[sc;table("optimized",in.drug.PS,'VariableNames',{'Source','PS_BBB_L_per_h'}); ...
 table("perfusion-limited",vs("perfusion_limited_PS_BBB"),'VariableNames',{'Source','PS_BBB_L_per_h'})];
end

function ind=buildIndividual(in,route)
phys=Physiology('human35m'); BW=double(getvalue(phys,'BW')); p=in.pk; d=in.drug; q=in.phys;
r=BW/p.BWref;
if route=="IV"
    absorption={'lambda_po',0/u.h,'F',0,'Tlag',0*u.h};
    dosing=Infusion('Paracetamol',0*u.h,in.doseIV*u.g,in.durationIV*u.min,'iv');
else
    absorption={'lambda_po',p.ka/u.h,'F',p.F,'Tlag',p.Tlag*u.h};
    dosing=Oral('Paracetamol',0*u.h,in.dosePO*u.g);
end
par=parameters(absorption{:},'V1',p.V1*r^p.eV1*u.L,'V2',p.V2*r^p.eV2*u.L, ...
 'V3',p.V3*r^p.eV3*u.L,'Q2',p.Q2*r^p.eQ2*u.L/u.h,'Q3',p.Q3*r^p.eQ3*u.L/u.h, ...
 'CL',p.CL*r^p.eCL*u.L/u.h,'fubrm',d.fuI,'Sex',1,'PS_BBB',d.PS*u.L/u.h, ...
 'V_brt',q.Vbrt*u.L,'V_brb',q.Vbrb*u.L,'endothelial_fraction',q.endo, ...
 'csf_fraction',q.csf,'cranial_csf_fraction',q.cranial,'brain_ecf_fraction',q.ecf, ...
 'Q_cbf',q.Qcbf*u.L/u.h,'Q_production',q.Qprod*u.L/u.h,'bulk_flow_fraction',q.bulk, ...
 'spinal_sink_fraction',q.sink,'spinal_outflow_fraction',q.out, ...
 'BCSFB_to_BBB_ratio',q.bcsfb,'BCM_to_BBB_ratio',q.bcm,'PS_e',d.PSE*u.L/u.h, ...
 'MW',d.MW*u.g/u.mol,'BP',d.BP,'fu_plasma',d.fuP,'fu_brain_plasma',d.fuBP, ...
 'fu_brain_ecf',d.fuE,'fu_brain_icf',d.fuI,'fu_csf',d.fuC, ...
 'fn_brain_plasma',d.fnBP,'fn_brain_ecf',d.fnE,'fn_brain_icf',d.fnI,'fn_csf',d.fnC);
ind=Individual('Virtual',1); ind.name='Human'; ind.physiology=phys; ind.dosing=dosing;
ind.drugdata=loaddrugdata('Paracetamol','species','human');
ind.sampling=Sampling((0:in.dt:in.tend)*u.h,FiveBrainobservables);
ind.model=paracetamol_Wang_plasma_5BRAIN; ind.model.par=par;
end

function out=simulateScenarios(ind,scenarios)
initialize(ind); simulate(ind); t=double(ind.observation.Time(1:8:end))/3600; n=numel(t);
v=double(ind.observation.Value)*1000; % model output g/L to mg/L
plasma=table(repmat("plasma",n,1),NaN(n,1),repmat("",n,1),t,v(1:8:end), ...
 'VariableNames',{'site','PS','source','time','conc'});
sites=["brain ECF","scsf"]; offsets=[4,8]; blocks=cell(height(scenarios)*2,1); k=0;
for s=1:2
 for i=1:height(scenarios)
  ind.model.par.PS_BBB=scenarios.PS_BBB_L_per_h(i)*u.L/u.h; initialize(ind); simulate(ind);
  values=double(ind.observation.Value)*1000; k=k+1;
  blocks{k}=table(repmat(sites(s),n,1),repmat(scenarios.PS_BBB_L_per_h(i),n,1), ...
   repmat(scenarios.Source(i),n,1),t,values(offsets(s):8:end), ...
   'VariableNames',{'site','PS','source','time','conc'});
 end
end
out=[plasma;vertcat(blocks{:})];
end

function value=exactValue(t,p)
row=t.Parameter==p; if ismember('ExactValue',t.Properties.VariableNames), col='ExactValue'; else, col='Value'; end
value=str2double(t.(col)(row));
end
function t=readStrings(file), o=detectImportOptions(file); o=setvartype(o,o.VariableNames,'string'); t=readtable(file,o); end
