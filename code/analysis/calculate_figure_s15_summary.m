%% Calculate the GMR and GMFE statistics used in Figure S15

clear; clc;
scriptDir=fileparts(mfilename('fullpath'));
repoRoot=fileparts(fileparts(scriptDir));
clinicalDir=fullfile(repoRoot,'data','clinical','figure_4');
simulationDir=fullfile(repoRoot,'results','simulations', ...
    'supplementary_figures','figure_s15');
outputDir=fullfile(repoRoot,'results','derived', ...
    'supplementary_figures','figure_s15');
if ~exist(outputDir,'dir'), mkdir(outputDir); end

compounds=["Ethanol","Diazepam","Paracetamol", ...
    "Ibuprofen","Indomethacin","Mannitol"];
blocks=cell(numel(compounds),1);
for c=1:numel(compounds)
 compound=compounds(c);
 clinicalFile=fullfile(clinicalDir,compound+"_clinical_data.csv");
 simulationFile=fullfile(simulationDir,compound+"_model_simulations.csv");
 assert(isfile(clinicalFile),'Missing clinical data: %s',clinicalFile);
 assert(isfile(simulationFile), ...
     'Missing simulations: %s. Run generate_figure_s15_model_simulations.m first.',simulationFile);
 clinical=readtable(clinicalFile,'TextType','string');
 simulations=readtable(simulationFile,'TextType','string');
 clinical.site(lower(clinical.site)=="scsf")="sCSF";
 clinical=clinical(clinical.site~="plasma",:);
 clinical.site(clinical.site=="brain mass")="brain ECF";
 simulations=simulations(simulations.site~="plasma",:);
 blocks{c}=calculateCompoundSummary(compound,clinical,simulations);
end
compoundSummary=vertcat(blocks{:});

sources=unique(compoundSummary.Source,'stable');
modelSummary=table('Size',[numel(sources),5], ...
    'VariableTypes',{'string','string','double','double','double'}, ...
    'VariableNames',{'Source','ModelClass','Mean_GMFE','SD_GMFE','N'});
for i=1:numel(sources)
 rows=compoundSummary.Source==sources(i);
 values=compoundSummary.GMFE(rows);
 modelSummary.Source(i)=sources(i);
 modelSummary.ModelClass(i)=compoundSummary.ModelClass(find(rows,1));
 modelSummary.Mean_GMFE(i)=mean(values,'omitnan');
 modelSummary.SD_GMFE(i)=std(values,'omitnan');
 modelSummary.N(i)=sum(isfinite(values));
end

compoundFile=fullfile(outputDir,'figure_s15_compound_summary.csv');
modelFile=fullfile(outputDir,'figure_s15_model_summary.csv');
writetable(compoundSummary,compoundFile); writetable(modelSummary,modelFile);
disp(compoundSummary); disp(modelSummary);
fprintf('Saved %s\nSaved %s\n',compoundFile,modelFile);

function output=calculateCompoundSummary(compound,clinical,simulations)
sources=unique(simulations.source,'stable'); sites=unique(simulations.site,'stable');
rows=cell(numel(sources)*numel(sites),1); k=0;
for i=1:numel(sources)
 source=string(sources(i)); sourceData=simulations(simulations.source==source,:);
 for j=1:numel(sites)
  site=string(sites(j)); sim=sourceData(sourceData.site==site,:);
  obs=clinical(clinical.site==site,:);
  if isempty(sim) || isempty(obs) || height(sim)<2, continue, end
  predicted=interp1(sim.time,sim.conc,obs.time,'linear','extrap');
  valid=predicted>0 & obs.conc>0 & isfinite(predicted) & isfinite(obs.conc);
  ratios=log(predicted(valid)./obs.conc(valid));
  if isempty(ratios), continue, end
  k=k+1;
  rows{k}=table(compound,site,source,string(sim.model_class(1)), ...
      string(sim.reference(1)),sim.PS(1),exp(mean(ratios)), ...
      exp(mean(abs(ratios))),sum(valid),'VariableNames', ...
      {'Compound','Site','Source','ModelClass','Reference','PS_BBB_L_per_h', ...
      'GMR','GMFE','N_observations'});
 end
end
rows=rows(1:k); output=vertcat(rows{:});
end
