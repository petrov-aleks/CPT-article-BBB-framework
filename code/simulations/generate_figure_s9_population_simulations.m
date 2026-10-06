%% Generate the population simulations used in Figure S9
% This workflow is computationally expensive. It intentionally does not run
% automatically from any figure-generation script.

clear; clc;
scriptDir = fileparts(mfilename('fullpath'));
repoRoot = fileparts(fileparts(scriptDir));
addpath(fullfile(scriptDir,'models'));

configureCompphysiol();
overwriteExisting = true; % Set true to regenerate every compound.
compounds = ["Ethanol","Diazepam","Paracetamol", ...
    "Ibuprofen","Indomethacin","Mannitol"];
outputDir = fullfile(repoRoot,'results','simulations', ...
    'supplementary_figures','figure_s9');
if ~exist(outputDir,'dir'), mkdir(outputDir); end

for compound = compounds
    outputFile = fullfile(outputDir,compound + "_population_simulations.csv");
    if isfile(outputFile) && ~overwriteExisting
        fprintf('Skipping %s: output already exists.\n',compound);
        continue
    end
    fprintf('Generating Figure S9 population simulation for %s...\n',compound);
    inputs = loadInputs(repoRoot,compound);
    rng(inputs.settings.RandomSeed,'twister');
    individuals = buildPopulation(inputs);
    initialize(individuals);
    applyVariableCNSFlows(individuals,inputs);
    simulate(individuals);
    summary = summarizePopulation(individuals,inputs.outputScale);
    writetable(summary,outputFile);
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

function inputs = loadInputs(repoRoot,compound)
files.s1 = fullfile(repoRoot,'data','model_parameters','table_s1_base_parameters.csv');
files.s1cv = fullfile(repoRoot,'data','model_parameters','table_s1_population_cv.csv');
files.s2 = fullfile(repoRoot,'data','model_parameters','table_s2_plasma_pk_parameters.csv');
files.s3 = fullfile(repoRoot,'results','tables',"table_s3_" + lower(compound) + ".csv");
files.s4 = fullfile(repoRoot,'data','clinical','table_s4_clinical_studies.csv');
files.settings = fullfile(repoRoot,'data','model_parameters','population_simulation_settings.csv');
for f = string(struct2cell(files))'
    if ~isfile(f), error('Required input was not found: %s',f); end
end
s1=readtable(files.s1,'TextType','string'); s1cv=readtable(files.s1cv,'TextType','string');
s2=readStrings(files.s2); s3=readStrings(files.s3);
s4=readtable(files.s4,'TextType','string'); settings=readtable(files.settings,'TextType','string');
v1=@(p) s1.Value(s1.Parameter==p); vcv=@(p) s1cv.Value(s1cv.Parameter==p);
ccv=@(p) s1cv.CV_percent(s1cv.Parameter==p);
valueColumn=compound+"_value"; cvColumn=compound+"_CV_percent";
v2=@(p) str2double(s2.(valueColumn)(s2.Parameter==p));
cv2=@(p) parseCV(s2.(cvColumn)(s2.Parameter==p));
v3=@(p) exactValue(s3,p);

inputs.compound=compound; inputs.s2value=v2; inputs.s2cv=cv2;
if any(compound==["Diazepam","Ibuprofen","Indomethacin"])
    inputs.outputScale=1e6; % g/L to ng/mL
elseif compound=="Paracetamol"
    inputs.outputScale=1e3; % g/L to mg/L
else
    inputs.outputScale=1;
end
inputs.physiology=struct('V_brt',v1("brain_tissue_volume"),'V_brb',v1("brain_blood_volume"), ...
 'endothelialFraction',v1("endothelial_fraction_of_total_brain"),'csfFraction',v1("csf_fraction_of_total_brain"), ...
 'cranialCSFFraction',v1("cranial_fraction_of_total_csf"),'brainECFFraction',v1("brain_ecf_fraction"), ...
 'Qcbf',v1("cerebral_blood_flow"),'Qproduction',v1("csf_production_rate"), ...
 'bulkFlowFraction',v1("bulk_flow_fraction"),'spinalSinkFraction',v1("spinal_sink_fraction"), ...
 'spinalOutflowFraction',v1("spinal_outflow_fraction"),'BCSFBtoBBBRatio',v1("bcsfb_to_bbb_surface_area_ratio"), ...
 'BCMtoBBBRatio',v1("bcm_to_bbb_surface_area_ratio"));
inputs.flow=struct('Qbulk',vcv("Q_bulk"),'QbulkCV',ccv("Q_bulk"), ...
 'Qssink',vcv("Q_ssink"),'QssinkCV',ccv("Q_ssink"), ...
 'Qsout',vcv("Q_sout"),'QsoutCV',ccv("Q_sout"));
inputs.drug=struct('MW',v3("MW"),'BP',v3("BP"),'fuPlasma',v3("fu_plasma"), ...
 'fuBrainPlasma',v3("fu_brain_plasma"),'fuBrainECF',v3("fu_brain_ECF"), ...
 'fuBrainICF',v3("fu_brain_ICF"),'fuCSF',v3("fu_CSF"), ...
 'fnBrainPlasma',v3("fn_brain_plasma"),'fnBrainECF',v3("fn_brain_ECF"), ...
 'fnBrainICF',v3("fn_brain_ICF"),'fnCSF',v3("fn_CSF"), ...
 'PSBBB',v3("PS_BBB"),'PSE',v3("PS_E"));
inputs.settings=settings(settings.Compound==compound,:);
rows=s4(s4.Drug==compound,:);
if compound=="Paracetamol", rows=rows(contains(rows.DosingRegimen,"IV"),:); end
inputs.doseValue=rows.ModelDoseValue(1); inputs.doseUnit=rows.ModelDoseUnit(1);
inputs.doseDuration=rows.ModelAdministrationDurationMin(1);
end

function individuals = buildPopulation(in)
N=in.settings.N; compound=in.compound; humanphys=Physiology('human35m');
BW=sampleTruncatedNormal(in.settings.BWMean_kg,in.settings.BWSD_kg, ...
    in.settings.BWMin_kg,in.settings.BWMax_kg,N);
pk=samplePK(in,BW); common=commonParameters(in);
individuals=Individual('Virtual',N); template=Individual('Virtual',1);
template.name='Human'; template.physiology=humanphys;
drugName=char(compound);
if compound=="Indomethacin", drugName='Indometacin'; end
template.drugdata=loaddrugdata(drugName,'species','human');
template.sampling=Sampling((0:in.settings.SamplingInterval_h:in.settings.SamplingEnd_h)*u.h,FiveBrainobservables);
switch compound
 case "Ethanol", template.model=ethanol_plasma_5BRAIN;
 case "Diazepam", template.model=diazepam_Hung_plasma_5BRAIN;
 case "Paracetamol", template.model=paracetamol_Wang_plasma_5BRAIN;
 case "Ibuprofen", template.model=ibuprofen_Gusthuys_plasma_5BRAIN;
 case "Indomethacin", template.model=indometacin_Saleh_plasma_5BRAIN;
 case "Mannitol", template.model=mannitol_plasma_5BRAIN;
end
for i=1:N
 individuals(i)=clone(template); individuals(i).dosing=makeDose(in,BW(i));
 switch compound
  case "Ethanol"
   sys={'V1',pk.V1(i)*u.L,'V2',pk.V2(i)*u.L,'Vmax',pk.Vmax(i)*u.g/u.h, ...
    'Km',pk.KM(i)*u.g/u.L,'Q',pk.Q(i)*u.L/u.h,'lambda_po',pk.ka(i)/u.h,'F',pk.F(i)};
  case "Diazepam"
   sys={'V1',pk.V1(i)*u.L,'V2',pk.V2(i)*u.L,'V3',pk.V3(i)*u.L, ...
    'CL1',pk.CL(i)*u.L/u.h,'CL2',pk.Q1(i)*u.L/u.h,'CL3',pk.Q2(i)*u.L/u.h, ...
    'lambda_po',pk.ka(i)/u.h,'F',pk.F(i),'Tlag',pk.Tlag(i)*u.h};
  case "Paracetamol"
   sys={'V1',pk.V1(i)*u.L,'V2',pk.V2(i)*u.L,'V3',pk.V3(i)*u.L, ...
    'Q2',pk.Q1(i)*u.L/u.h,'Q3',pk.Q2(i)*u.L/u.h,'CL',pk.CL(i)*u.L/u.h, ...
    'lambda_po',0/u.h,'F',0,'Tlag',0*u.h,'fubrm',in.drug.fuBrainICF,'Sex',1};
  case "Ibuprofen"
   sys={'V1',pk.V1(i)*u.L,'CL',pk.CL(i)*u.L/u.h,'lambda_po',pk.ka(i)/u.h, ...
    'F',pk.F(i),'Tlag',pk.Tlag(i)*u.h};
  case "Indomethacin"
   sys={'V1',pk.V1(i)*u.L,'V2',pk.V2(i)*u.L,'CL1',pk.CL(i)*u.L/u.h, ...
    'CL2',pk.Q1(i)*u.L/u.h,'lambda_po',pk.ka(i)/u.h,'F',pk.F(i),'Tlag',pk.Tlag(i)*u.h};
  case "Mannitol"
   sys={'V1',pk.V1(i)*u.L,'V2',pk.V2(i)*u.L,'Q',pk.Q(i)*u.L/u.h,'CL',pk.CL(i)*u.L/u.h, ...
    'lambda_po',0/u.h,'F',0,'Tlag',0*u.h};
 end
 individuals(i).model.par=parameters(sys{:},common{:});
end
end

function pk=samplePK(in,BW)
N=numel(BW); v=in.s2value; cv=in.s2cv; c=in.compound;
s=@(p) sampleLognormal(v(p),cv(p),N); fixed=@(p) repmat(v(p),1,N);
switch c
 case "Ethanol"
  pk.V1=s("Vcntr"); pk.V2=s("Vpone").*(BW/v("BWref")).^v("BWVpone");
  pk.Q=s("Qcpone"); pk.KM=fixed("KM"); pk.ka=s("ka"); pk.F=fixed("F");
  food=s("FEVmax"); iov=sampleLognormal(1,in.settings.VmaxIOV_CV_percent,N);
  pk.Vmax=s("Vmax").*food.*iov;
 case "Diazepam"
  pk.V1=s("Vcntr"); pk.V2=s("Vpone"); pk.V3=s("Vptwo"); pk.CL=s("CL");
  pk.Q1=s("Qcpone"); pk.Q2=s("Qcptwo"); pk.ka=s("ka"); pk.F=fixed("F"); pk.Tlag=fixed("Tlag");
 case "Paracetamol"
  ref=v("BWref"); pk.V1=s("Vcntr").*(BW/ref).^v("BWVcntr");
  pk.V2=s("Vpone").*(BW/ref).^v("BWVpone"); pk.V3=s("Vptwo").*(BW/ref).^v("BWVptwo");
  pk.Q1=s("Qcpone").*(BW/ref).^v("BWQcpone"); pk.Q2=s("Qcptwo").*(BW/ref).^v("BWQcptwo");
  pk.CL=s("CL").*(BW/ref).^v("BWCL");
 case "Ibuprofen"
  ref=v("BWref"); pk.V1=s("Vcntr").*(BW/ref).^v("BWVcntr");
  pk.CL=s("CL").*(BW/ref).^v("BWCL"); pk.ka=s("ka"); pk.F=fixed("F"); pk.Tlag=s("Tlag");
 case "Indomethacin"
  pk.V1=s("Vcntr"); pk.V2=s("Vpone"); pk.CL=s("CL"); pk.Q1=s("Qcpone");
  pk.ka=s("ka"); pk.F=fixed("F"); pk.Tlag=fixed("Tlag");
 case "Mannitol"
  names=["Vcntr","Vpone","Qcpone","CL"];
  means=arrayfun(@(p)v(p),names); cvs=arrayfun(@(p)cv(p),names)/100;
  vars=log(1+cvs.^2); a=vars(1); b=vars(3)-a; cc=vars(2)-vars(3); d=vars(4)-a;
  Sigma=[a,a,a,a; a,a+b,a+b,a; a,a+b,a+b+cc,a; a,a,a,a+d];
  mu=log(means)-0.5*diag(Sigma)'; Z=randn(4,N); X=exp(mu'+chol(Sigma,'lower')*Z);
  pk.V1=X(1,:); pk.V2=X(2,:); pk.Q=X(3,:); pk.CL=X(4,:);
end
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
c=in.compound; value=in.doseValue; duration=in.doseDuration*u.min;
if c=="Ethanol"
    dose=Infusion('Ethanol',0*u.h,value*u.g/u.kg*bw*u.kg,duration,'drink');
elseif c=="Paracetamol"
    dose=Infusion('Paracetamol',0*u.h,value*u.g,duration,'iv');
elseif c=="Mannitol"
    dose=Infusion('Mannitol',0*u.h,value*u.mg/u.kg*bw*u.kg,duration,'iv');
elseif c=="Diazepam"
    dose=Oral('Diazepam',0*u.h,value*u.mg);
elseif c=="Ibuprofen"
    dose=Oral('Ibuprofen',0*u.h,value*u.mg);
else
    dose=Oral('Indometacin',0*u.h,value*u.mg);
end
end

function applyVariableCNSFlows(individuals,in)
N=numel(individuals); q=in.flow;
bulk=sampleLognormal(q.Qbulk,q.QbulkCV,N); sink=sampleLognormal(q.Qssink,q.QssinkCV,N);
out=sampleLognormal(q.Qsout,q.QsoutCV,N); production=in.physiology.Qproduction;
for i=1:N
 individuals(i).model.setup.Q.bulk=bulk(i)*u.L/u.h;
 individuals(i).model.setup.Q.ssink=sink(i)*u.L/u.h;
 individuals(i).model.setup.Q.sout=out(i)*u.L/u.h;
 individuals(i).model.setup.Q.csink=(production-sink(i))*u.L/u.h;
 individuals(i).model.setup.Q.sin=(sink(i)+out(i))*u.L/u.h;
end
end

function out=summarizePopulation(individuals,outputScale)
sites=["plasma","brint","sCSF"]; offsets=[1,4,8]; N=numel(individuals); blocks=cell(3,1);
for s=1:3
 t=double(individuals(1).observation.Time(offsets(s):8:end))/3600; values=zeros(numel(t),N);
 for i=1:N
  values(:,i)=double(individuals(i).observation.Value(offsets(s):8:end))*outputScale;
 end
 if s>1 && ~strcmp(sites(s),"brint"), values=values; end %#ok<ASGSL>
 pct=prctile(values,[5 50 95],2);
 blocks{s}=table(repmat(sites(s),numel(t),1),t,pct(:,1),pct(:,2),pct(:,3), ...
  'VariableNames',{'sim','time','p5','p50','p95'});
end
out=vertcat(blocks{:});
end

function x=sampleLognormal(typical,cvPercent,N)
if isnan(cvPercent)||cvPercent==0, x=repmat(typical,1,N); return; end
omega=sqrt(log(1+(cvPercent/100)^2)); x=exp(log(typical)+omega*randn(1,N));
end
function x=sampleTruncatedNormal(mu,sd,lo,hi,N)
if sd==0, x=repmat(mu,1,N); else, x=random(truncate(makedist('Normal','mu',mu,'sigma',sd),lo,hi),1,N); end
end
function value=parseCV(x), value=str2double(x); if isnan(value), value=0; end, end
function value=exactValue(t,p)
row=t.Parameter==p; if ismember('ExactValue',t.Properties.VariableNames), col='ExactValue'; else, col='Value'; end
value=str2double(t.(col)(row));
end
function t=readStrings(file)
o=detectImportOptions(file); o=setvartype(o,o.VariableNames,'string'); t=readtable(file,o);
end
