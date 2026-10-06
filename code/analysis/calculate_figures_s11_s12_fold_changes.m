%% Calculate ratios and fold errors relative to optimized simulations

clear; clc;
scriptDir=fileparts(mfilename('fullpath'));
repoRoot=fileparts(fileparts(scriptDir));
inputFile=fullfile(repoRoot,'results','summary_data','figures_s11_s12_exposure_metrics.csv');
if ~isfile(inputFile)
    error('Run calculate_figures_s11_s12_exposure_metrics.m first.');
end
metrics=readtable(inputFile,'TextType','string');
groups=unique(metrics(:,{'Compound','Site'}),'rows','stable');
blocks=cell(0,1);

for g=1:height(groups)
    rows=metrics(metrics.Compound==groups.Compound(g) & metrics.Site==groups.Site(g),:);
    reference=rows(rows.Source=="optimized",:);
    assert(height(reference)==1,'Expected one optimized scenario for %s, %s.', ...
        groups.Compound(g),groups.Site(g));
    ratioCmax=rows.Cmax/reference.Cmax;
    ratioAUC=rows.AUC/reference.AUC;
    foldCmax=max(ratioCmax,1./ratioCmax);
    foldAUC=max(ratioAUC,1./ratioAUC);
    n=height(rows);
    blocks{end+1,1}=table([rows.Compound;rows.Compound],[rows.Site;rows.Site], ...
        [rows.Source;rows.Source],[repmat("Cmax",n,1);repmat("AUC",n,1)], ...
        [ratioCmax;ratioAUC],[foldCmax;foldAUC], ...
        'VariableNames',{'Compound','Site','Source','Metric','Ratio','FoldError'}); %#ok<SAGROW>
end

foldChanges=vertcat(blocks{:});
foldChanges=sortrows(foldChanges,{'Compound','Site','Metric','Source'});
outputFile=fullfile(repoRoot,'results','summary_data','figures_s11_s12_fold_changes.csv');
writetable(foldChanges,outputFile);
fprintf('Fold-change summary saved: %s\n',outputFile);
