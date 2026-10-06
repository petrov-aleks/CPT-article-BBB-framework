%% Calculate Cmax and AUC summary data for Figures S11 and S12

clear; clc;
scriptDir=fileparts(mfilename('fullpath'));
repoRoot=fileparts(fileparts(scriptDir));
generatedDir=fullfile(repoRoot,'results','simulations','main_figures','figure_4');
outputDir=fullfile(repoRoot,'results','summary_data');
if ~exist(outputDir,'dir'), mkdir(outputDir); end

compounds=["Ethanol","Diazepam","Paracetamol","Ibuprofen","Indomethacin","Mannitol"];
sites=["brain ECF","scsf"];
blocks=cell(0,1);

for compound=compounds
    file=compound+"_simulations.csv";
    inputFile=fullfile(generatedDir,file);
    if ~isfile(inputFile)
        error('Missing Figure 4 simulations for %s. Run its simulation workflow first.',compound);
    end
    opts=detectImportOptions(inputFile); opts=setvartype(opts,{'site','source'},'string');
    sim=readtable(inputFile,opts);
    sim.source(sim.source=="in vitro MDCK/UC" | sim.source=="in vitro Caco-2")="in vitro cell-based";
    for site=sites
        siteData=sim(strcmpi(sim.site,site) & strlength(sim.source)>0,:);
        sources=unique(siteData.source,'stable');
        for source=sources'
            rows=siteData(siteData.source==source,:);
            rows=sortrows(rows,'time');
            cmax=max(rows.conc);
            auc=trapz(rows.time,rows.conc);
            blocks{end+1,1}=table(compound,site,source,cmax,auc, ...
                'VariableNames',{'Compound','Site','Source','Cmax','AUC'}); %#ok<SAGROW>
        end
    end
end

metrics=vertcat(blocks{:});
outputFile=fullfile(outputDir,'figures_s11_s12_exposure_metrics.csv');
writetable(metrics,outputFile);
fprintf('Exposure metrics saved: %s\n',outputFile);
