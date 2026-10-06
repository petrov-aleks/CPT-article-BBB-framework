%% Generate the individual surface-area simulations used in Figure S10

clear; clc;
scriptDir = fileparts(mfilename('fullpath'));
repoRoot = fileparts(fileparts(scriptDir));
addpath(fullfile(scriptDir,'models'));
configureCompphysiol();

compounds = ["Ethanol","Diazepam","Paracetamol", ...
    "Ibuprofen","Indomethacin","Mannitol"];
scenarios = readtable(fullfile(repoRoot,'data','model_parameters', ...
    'figure_s10_surface_area_scenarios.csv'),'TextType','string');
outputDir = fullfile(repoRoot,'results','simulations', ...
    'supplementary_figures','figure_s10');
if ~exist(outputDir,'dir'), mkdir(outputDir); end

for compound = compounds
    fprintf('Generating Figure S10 simulations for %s...\n',compound);
    inputs = loadInputs(repoRoot,compound);
    individual = buildIndividual(inputs);
    simulations = simulateSurfaceAreas(individual,scenarios,inputs.drug.PSBBB,inputs.outputScale);
    outputFile = fullfile(outputDir,compound + "_surface_area_simulations.csv");
    writetable(simulations,outputFile);
    fprintf('Saved %s\n',outputFile);
end

%% Local functions
function configureCompphysiol()
setoptcompphysiol('DrugDB',@initdrugdbCNS);
setoptcompphysiol('DrugTemplate',@drugtemplateCNS);
setoptcompphysiol('ObservableTemplate',@obstemplateCNS);
setoptcompphysiol('PlotTemplate',@plottemplateCNS);
setoptcompphysiol('DisplayUnits',{'g/L','L/h'});
end

function in = loadInputs(repoRoot,compound)
files.s1=fullfile(repoRoot,'data','model_parameters','table_s1_base_parameters.csv');
files.s2=fullfile(repoRoot,'data','model_parameters','table_s2_plasma_pk_parameters.csv');
files.s3=fullfile(repoRoot,'results','tables',"table_s3_"+lower(compound)+".csv");
files.s4=fullfile(repoRoot,'data','clinical','table_s4_clinical_studies.csv');
files.settings=fullfile(repoRoot,'data','model_parameters',lower(compound)+"_simulation_settings.csv");
for f=string(struct2cell(files))', if ~isfile(f), error('Required input was not found: %s',f); end, end
s1=readtable(files.s1,'TextType','string'); s2=readStrings(files.s2);
s3=readStrings(files.s3); s4=readtable(files.s4,'TextType','string'); settings=readStrings(files.settings);
v1=@(p)s1.Value(s1.Parameter==p); vc=compound+"_value"; cc=compound+"_CV_percent";
v2=@(p)str2double(s2.(vc)(s2.Parameter==p)); cv2=@(p)parseCV(s2.(cc)(s2.Parameter==p));
v3=@(p)exactValue(s3,p); vs=@(p)str2double(settings.Value(settings.Parameter==p));
in.compound=compound; in.s2value=v2; in.s2cv=cv2;
if any(compound==["Diazepam","Ibuprofen","Indomethacin"])
    in.outputScale=1e6; % g/L to ng/mL
elseif compound=="Paracetamol"
    in.outputScale=1e3; % g/L to mg/L
else
    in.outputScale=1;
end
in.physiology=struct('V_brt',v1("brain_tissue_volume"),'V_brb',v1("brain_blood_volume"), ...
 'endothelialFraction',v1("endothelial_fraction_of_total_brain"),'csfFraction',v1("csf_fraction_of_total_brain"), ...
 'cranialCSFFraction',v1("cranial_fraction_of_total_csf"),'brainECFFraction',v1("brain_ecf_fraction"), ...
 'Qcbf',v1("cerebral_blood_flow"),'Qproduction',v1("csf_production_rate"), ...
 'bulkFlowFraction',v1("bulk_flow_fraction"),'spinalSinkFraction',v1("spinal_sink_fraction"), ...
 'spinalOutflowFraction',v1("spinal_outflow_fraction"),'BCSFBtoBBBRatio',v1("bcsfb_to_bbb_surface_area_ratio"), ...
 'BCMtoBBBRatio',v1("bcm_to_bbb_surface_area_ratio"));
in.drug=struct('MW',v3("MW"),'BP',v3("BP"),'fuPlasma',v3("fu_plasma"), ...
 'fuBrainPlasma',v3("fu_brain_plasma"),'fuBrainECF',v3("fu_brain_ECF"), ...
 'fuBrainICF',v3("fu_brain_ICF"),'fuCSF',v3("fu_CSF"), ...
 'fnBrainPlasma',v3("fn_brain_plasma"),'fnBrainECF',v3("fn_brain_ECF"), ...
 'fnBrainICF',v3("fn_brain_ICF"),'fnCSF',v3("fn_CSF"), ...
 'PSBBB',v3("PS_BBB"),'PSE',v3("PS_E"));
rows=s4(s4.Drug==compound,:); if any(compound==["Paracetamol","Mannitol"]), rows=rows(contains(rows.DosingRegimen,"IV"),:); end
in.doseValue=rows.ModelDoseValue(1); in.doseDuration=rows.ModelAdministrationDurationMin(1);
in.samplingInterval=vs("sampling_interval");
ends=containers.Map({'Ethanol','Diazepam','Paracetamol','Ibuprofen','Indomethacin','Mannitol'},[10,25,14,25,25,7]);
in.samplingEnd=ends(char(compound));
end

function ind=buildIndividual(in)
c=in.compound; v=in.s2value; cv=in.s2cv; phys=Physiology('human35m'); BW=double(getvalue(phys,'BW'));
common=commonParameters(in); typical=@(p)v(p)/sqrt(1+(cv(p)/100)^2);
switch c
 case "Ethanol"
  sys={'V1',v("Vcntr")*u.L,'V2',v("Vpone")*(BW/v("BWref"))^v("BWVpone")*u.L, ...
   'Vmax',v("Vmax")*v("FEVmax")*u.g/u.h,'Km',v("KM")*u.g/u.L,'Q',v("Qcpone")*u.L/u.h, ...
   'lambda_po',v("ka")/u.h,'F',v("F")}; model=ethanol_plasma_5BRAIN;
 case "Diazepam"
  sys={'V1',v("Vcntr")*u.L,'V2',v("Vpone")*u.L,'V3',v("Vptwo")*u.L, ...
   'CL1',v("CL")*u.L/u.h,'CL2',v("Qcpone")*u.L/u.h,'CL3',v("Qcptwo")*u.L/u.h, ...
   'lambda_po',v("ka")/u.h,'F',v("F"),'Tlag',v("Tlag")*u.h}; model=diazepam_Hung_plasma_5BRAIN;
 case "Paracetamol"
  r=v("BWref"); sys={'V1',v("Vcntr")*(BW/r)^v("BWVcntr")*u.L, ...
   'V2',v("Vpone")*(BW/r)^v("BWVpone")*u.L,'V3',v("Vptwo")*(BW/r)^v("BWVptwo")*u.L, ...
   'Q2',v("Qcpone")*(BW/r)^v("BWQcpone")*u.L/u.h,'Q3',v("Qcptwo")*(BW/r)^v("BWQcptwo")*u.L/u.h, ...
   'CL',v("CL")*(BW/r)^v("BWCL")*u.L/u.h,'lambda_po',0/u.h,'F',0,'Tlag',0*u.h, ...
   'fubrm',in.drug.fuBrainICF,'Sex',1}; model=paracetamol_Wang_plasma_5BRAIN;
 case "Ibuprofen"
  r=v("BWref"); sys={'V1',v("Vcntr")*(BW/r)^v("BWVcntr")*u.L, ...
   'CL',v("CL")*(BW/r)^v("BWCL")*u.L/u.h,'lambda_po',v("ka")/u.h,'F',v("F"),'Tlag',v("Tlag")*u.h};
  model=ibuprofen_Gusthuys_plasma_5BRAIN;
 case "Indomethacin"
  sys={'V1',v("Vcntr")*u.L,'V2',v("Vpone")*u.L,'CL1',v("CL")*u.L/u.h, ...
   'CL2',v("Qcpone")*u.L/u.h,'lambda_po',v("ka")/u.h,'F',v("F"),'Tlag',v("Tlag")*u.h};
  model=indometacin_Saleh_plasma_5BRAIN;
 case "Mannitol"
  sys={'V1',typical("Vcntr")*u.L,'V2',typical("Vpone")*u.L, ...
   'Q',typical("Qcpone")*u.L/u.h,'CL',typical("CL")*u.L/u.h, ...
   'lambda_po',0/u.h,'F',0,'Tlag',0*u.h}; model=mannitol_plasma_5BRAIN;
end
ind=Individual('Virtual',1); ind.name='Human'; ind.physiology=phys; ind.dosing=makeDose(in,BW);
drugName=char(c); if c=="Indomethacin", drugName='Indometacin'; end
ind.drugdata=loaddrugdata(drugName,'species','human');
ind.sampling=Sampling((0:in.samplingInterval:in.samplingEnd)*u.h,FiveBrainobservables);
ind.model=model; ind.model.par=parameters(sys{:},common{:});
end

function pairs=commonParameters(in)
p=in.physiology; d=in.drug;
pairs={'PS_BBB',d.PSBBB*u.L/u.h,'V_brt',p.V_brt*u.L,'V_brb',p.V_brb*u.L, ...
 'endothelial_fraction',p.endothelialFraction,'csf_fraction',p.csfFraction, ...
 'cranial_csf_fraction',p.cranialCSFFraction,'brain_ecf_fraction',p.brainECFFraction, ...
 'Q_cbf',p.Qcbf*u.L/u.h,'Q_production',p.Qproduction*u.L/u.h, ...
 'bulk_flow_fraction',p.bulkFlowFraction,'spinal_sink_fraction',p.spinalSinkFraction, ...
 'spinal_outflow_fraction',p.spinalOutflowFraction,'BCSFB_to_BBB_ratio',p.BCSFBtoBBBRatio, ...
 'BCM_to_BBB_ratio',p.BCMtoBBBRatio,'PS_e',d.PSE*u.L/u.h,'MW',d.MW*u.g/u.mol, ...
 'BP',d.BP,'fu_plasma',d.fuPlasma,'fu_brain_plasma',d.fuBrainPlasma, ...
 'fu_brain_ecf',d.fuBrainECF,'fu_brain_icf',d.fuBrainICF,'fu_csf',d.fuCSF, ...
 'fn_brain_plasma',d.fnBrainPlasma,'fn_brain_ecf',d.fnBrainECF, ...
 'fn_brain_icf',d.fnBrainICF,'fn_csf',d.fnCSF};
end

function dose=makeDose(in,bw)
c=in.compound; x=in.doseValue; duration=in.doseDuration*u.min;
if c=="Ethanol"
    dose=Infusion('Ethanol',0*u.h,x*u.g/u.kg*bw*u.kg,duration,'drink');
elseif c=="Paracetamol"
    dose=Infusion('Paracetamol',0*u.h,x*u.g,duration,'iv');
elseif c=="Mannitol"
    dose=Infusion('Mannitol',0*u.h,x*u.mg/u.kg*bw*u.kg,duration,'iv');
elseif c=="Diazepam"
    dose=Oral('Diazepam',0*u.h,x*u.mg);
elseif c=="Ibuprofen"
    dose=Oral('Ibuprofen',0*u.h,x*u.mg);
else
    dose=Oral('Indometacin',0*u.h,x*u.mg);
end
end

function out=simulateSurfaceAreas(ind,scenarios,optimizedPS,outputScale)
% The optimized value refers to 15 m2 BBB area; preserve intrinsic P0.
p0=optimizedPS/15; sites=["brain ECF","brain ICF","sCSF"]; offsets=[4,5,8];
blocks=cell(height(scenarios)*3,1); k=0;
for i=1:height(scenarios)
 initialize(ind);
 ind.model.setup.PS.b=p0*scenarios.SA_BBB_m2(i)*u.L/u.h;
 ind.model.setup.PS.c=p0*scenarios.SA_BCSFB_m2(i)*u.L/u.h;
 ind.model.setup.PS.i=p0*scenarios.SA_BCM_m2(i)*u.L/u.h;
 simulate(ind); t=double(ind.observation.Time(1:8:end))/3600;
 values=double(ind.observation.Value)*outputScale;
 for s=1:3
  k=k+1; blocks{k}=table(repmat(sites(s),numel(t),1),t,values(offsets(s):8:end), ...
   repmat(scenarios.Scenario(i),numel(t),1),'VariableNames',{'site','time','conc','PS'});
 end
end
out=vertcat(blocks{:});
end

function value=parseCV(x), value=str2double(x); if isnan(value), value=0; end, end
function value=exactValue(t,p)
row=t.Parameter==p; if ismember('ExactValue',t.Properties.VariableNames), col='ExactValue'; else, col='Value'; end
value=str2double(t.(col)(row));
end
function t=readStrings(file)
o=detectImportOptions(file); o=setvartype(o,o.VariableNames,'string'); t=readtable(file,o);
end
