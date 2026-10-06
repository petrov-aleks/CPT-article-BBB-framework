%% Generate Figure S14: comparison of alternative in silico models

clear; clc; close all;

scriptDir = fileparts(mfilename('fullpath'));
repoRoot = fileparts(fileparts(scriptDir));
modelFile = fullfile(repoRoot,'results','tables','table_s5.csv');
trainingFile = fullfile(repoRoot,'data','model_development', ...
    'in_silico_model_training_data.csv');
outputDir = fullfile(repoRoot,'results','figures','supplementary');
summaryDir = fullfile(repoRoot,'results','summary_data');

assert(isfile(modelFile), ...
    'Table S5 was not found. Run code/tables/generate_table_s5.m first.');
assert(isfile(trainingFile),'Training data were not found: %s',trainingFile);
if ~exist(outputDir,'dir'), mkdir(outputDir); end
if ~exist(summaryDir,'dir'), mkdir(summaryDir); end

models = readtable(modelFile,'TextType','string');
training = readtable(trainingFile,'TextType','string');
requiredModelColumns = ["Model","Intercept","LogP_Coefficient", ...
    "Log10_MW_Coefficient","MW_over_100_Coefficient"];
assert(all(ismember(requiredModelColumns,string(models.Properties.VariableNames))), ...
    'Table S5 lacks the machine-readable model coefficients.');

smith = selectModel(models,"Smith et al. 2024");
luco = selectModel(models,"Luco et al. 2006");
tsinman = selectModel(models,"Tsinman-based model");
grumetto = selectModel(models,"Grumetto et al. 2016");

% Figure settings retained from the manuscript workflow.
figureWidth = 17.8; figureHeight = 15;
axisLineWidth = 1.2; modelLineWidth = 1.5;
grumettoLineWidth = 1.8; dataLineWidth = 1.2; markerSize = 5;
tickFont = 10; labelFont = 10; titleFont = 10; legendFont = 10;
mwLines = [50 200 400];
lineStyles = {'--','-',':'};
logPLimits = [-6 8]; logPTicks = -6:2:8;
colors.Smith = [0.9255 0.0039 0.0039];
colors.Luco = [0.6510 0.0039 0.4745];
colors.Tsinman = [0 0.4745 0];
colors.Grumetto = [0.9961 0.7059 0.3843];

brainWeight_g = 1450;
bbbSurfaceArea_m2 = 15;
ps0ToPSBBB = brainWeight_g*3600/1000;
p0ToPSBBB = 0.01*bbbSurfaceArea_m2*3600*1000;
logPGrid = linspace(logPLimits(1),logPLimits(2),300).';

fig = figure('Color','w','Units','centimeters', ...
    'Position',[2 2 figureWidth figureHeight]);
layout = tiledlayout(fig,2,2,'TileSpacing','compact','Padding','compact');

% A. Perfusion models and their training datasets.
ax1 = nexttile(layout,1); hold(ax1,'on'); box(ax1,'on');
plotTraining(ax1,training,"Smith_2024",colors.Smith,markerSize,dataLineWidth);
plotTraining(ax1,training,"Luco_2006",colors.Luco,markerSize,dataLineWidth);
for i = 1:numel(mwLines)
    plot(ax1,logPGrid,evaluateLogModel(smith,logPGrid,mwLines(i)), ...
        'Color',colors.Smith,'LineStyle',lineStyles{i}, ...
        'LineWidth',modelLineWidth,'HandleVisibility','off');
    plot(ax1,logPGrid,evaluateLogModel(luco,logPGrid,mwLines(i)), ...
        'Color',colors.Luco,'LineStyle',lineStyles{i}, ...
        'LineWidth',modelLineWidth,'HandleVisibility','off');
end
formatAxis(ax1,tickFont,axisLineWidth,logPLimits,logPTicks);
ylim(ax1,[-8 8]); yticks(ax1,-8:2:8);
xlabel(ax1,'logP','FontSize',labelFont);
ylabel(ax1,'logPS_0','FontSize',labelFont);
title(ax1,'A. Smith 2024 vs Luco 2006','FontSize',titleFont,'FontWeight','bold');

% B. PAMPA models and their training datasets.
ax2 = nexttile(layout,2); hold(ax2,'on'); box(ax2,'on');
plotTraining(ax2,training,"Tsinman_2011",colors.Tsinman,markerSize,dataLineWidth);
plotTraining(ax2,training,"Grumetto_2016",colors.Grumetto,markerSize,dataLineWidth);
for i = 1:numel(mwLines)
    plot(ax2,logPGrid,evaluateLogModel(tsinman,logPGrid,mwLines(i)), ...
        'Color',colors.Tsinman,'LineStyle',lineStyles{i}, ...
        'LineWidth',modelLineWidth,'HandleVisibility','off');
end
plot(ax2,logPGrid,evaluateLogModel(grumetto,logPGrid,NaN), ...
    'Color',colors.Grumetto,'LineStyle','-', ...
    'LineWidth',grumettoLineWidth,'HandleVisibility','off');
formatAxis(ax2,tickFont,axisLineWidth,logPLimits,logPTicks);
ylim(ax2,[-10 2]); yticks(ax2,-10:2:2);
xlabel(ax2,'logP','FontSize',labelFont);
ylabel(ax2,'logP_0','FontSize',labelFont);
title(ax2,'B. Tsinman 2011 vs Grumetto 2016', ...
    'FontSize',titleFont,'FontWeight','bold');

% C. All equations translated to the common PS_BBB scale.
ax3 = nexttile(layout,3); hold(ax3,'on'); box(ax3,'on');
summary = table();
for i = 1:numel(mwLines)
    mw = mwLines(i);
    summary = [summary; plotPSBBB(ax3,smith,"Smith",logPGrid,mw, ...
        ps0ToPSBBB,colors.Smith,lineStyles{i},modelLineWidth)]; %#ok<AGROW>
    summary = [summary; plotPSBBB(ax3,luco,"Luco",logPGrid,mw, ...
        ps0ToPSBBB,colors.Luco,lineStyles{i},modelLineWidth)]; %#ok<AGROW>
    summary = [summary; plotPSBBB(ax3,tsinman,"Tsinman",logPGrid,mw, ...
        p0ToPSBBB,colors.Tsinman,lineStyles{i},modelLineWidth)]; %#ok<AGROW>
end
summary = [summary; plotPSBBB(ax3,grumetto,"Grumetto",logPGrid,NaN, ...
    p0ToPSBBB,colors.Grumetto,'-',grumettoLineWidth)];
set(ax3,'YScale','log');
formatAxis(ax3,tickFont,axisLineWidth,logPLimits,logPTicks);
ylim(ax3,[1e-6 1e10]); yticks(ax3,10.^(-6:2:10));
xlabel(ax3,'logP','FontSize',labelFont);
ylabel(ax3,'PS_{BBB} (L/h)','FontSize',labelFont);
title(ax3,'C. All model predictions translated to PS_{BBB}', ...
    'FontSize',titleFont,'FontWeight','bold');

% Permeability-limited region used in the manuscript.
permeabilityBoundary = 5; zoneColor = [1 0.5 0];
zone = patch(ax3,[logPLimits(1) logPLimits(2) logPLimits(2) logPLimits(1)], ...
    [1e-6 1e-6 permeabilityBoundary permeabilityBoundary],zoneColor, ...
    'FaceAlpha',0.08,'EdgeColor','none','HandleVisibility','off');
uistack(zone,'bottom');
yline(ax3,permeabilityBoundary,'--','Color',zoneColor, ...
    'LineWidth',0.8,'HandleVisibility','off');
text(ax3,logPLimits(2)-0.2,10^(-6+0.08*16), ...
    'Permeability-limited (E% < 10)','Color',zoneColor, ...
    'HorizontalAlignment','right','FontWeight','bold','FontSize',10);

% Legends occupy the fourth tile, as in the manuscript figure.
legendTile = nexttile(layout,4); axis(legendTile,'off'); drawnow;
position = legendTile.Position; delete(legendTile);
modelPosition = [position(1),position(2),position(3)*0.62,position(4)];
dataPosition = [position(1)+position(3)*0.64,position(2),position(3)*0.36,position(4)];

modelAxis = axes(fig,'Position',modelPosition,'Visible','off'); hold(modelAxis,'on');
[modelHandles,modelLabels] = modelLegend(modelAxis,mwLines,lineStyles,colors, ...
    modelLineWidth,grumettoLineWidth);
modelLegendObject = legend(modelAxis,modelHandles,modelLabels, ...
    'Location','northwest','Box','on','FontSize',legendFont);
modelLegendObject.Title.String = 'Models';

dataAxis = axes(fig,'Position',dataPosition,'Visible','off'); hold(dataAxis,'on');
names = ["Smith","Luco","Tsinman","Grumetto"];
dataHandles = gobjects(numel(names),1);
for i = 1:numel(names)
    dataHandles(i) = plot(dataAxis,nan,nan,'o','MarkerSize',markerSize, ...
        'MarkerFaceColor','none','MarkerEdgeColor',colors.(names(i)), ...
        'LineStyle','none','LineWidth',dataLineWidth);
end
dataLegendObject = legend(dataAxis,dataHandles,cellstr(names), ...
    'Location','northwest','Box','on','FontSize',legendFont);
dataLegendObject.Title.String = 'Data';

pdfFile = fullfile(outputDir,'Figure_S14.pdf');
pngFile = fullfile(outputDir,'Figure_S14.png');
summaryFile = fullfile(summaryDir,'figure_s14_model_trajectories.csv');
exportgraphics(fig,pdfFile,'ContentType','vector');
exportgraphics(fig,pngFile,'Resolution',600);
writetable(summary,summaryFile);
fprintf('Saved %s\nSaved %s\nSaved %s\n',pdfFile,pngFile,summaryFile);

function row = selectModel(models,name)
row = models(models.Model==name,:);
assert(height(row)==1,'Expected exactly one Table S5 row for %s.',name);
end

function values = evaluateLogModel(model,logP,mw)
values = model.Intercept + model.LogP_Coefficient.*logP;
if model.Log10_MW_Coefficient ~= 0
    assert(isfinite(mw) && mw>0,'A positive molecular weight is required.');
    values = values + model.Log10_MW_Coefficient.*log10(mw);
end
if model.MW_over_100_Coefficient ~= 0
    assert(isfinite(mw) && mw>0,'A positive molecular weight is required.');
    values = values + model.MW_over_100_Coefficient.*(mw/100);
end
end

function plotTraining(ax,data,source,color,markerSize,lineWidth)
subset = data(data.source_model==source,:);
assert(~isempty(subset),'No training data found for %s.',source);
plot(ax,subset.logP,subset.observed_log_value,'o','MarkerSize',markerSize, ...
    'MarkerFaceColor','none','MarkerEdgeColor',color,'LineStyle','none', ...
    'LineWidth',lineWidth,'HandleVisibility','off');
end

function formatAxis(ax,fontSize,lineWidth,xLimits,xTicks)
ax.LineWidth = lineWidth; ax.XAxis.FontSize = fontSize; ax.YAxis.FontSize = fontSize;
xlim(ax,xLimits); xticks(ax,xTicks);
end

function output = plotPSBBB(ax,model,name,logP,mw,conversion,color,style,width)
logPrediction = evaluateLogModel(model,logP,mw);
psbbb = 10.^logPrediction.*conversion;
plot(ax,logP,psbbb,'Color',color,'LineStyle',style, ...
    'LineWidth',width,'HandleVisibility','off');
output = table(repmat(name,numel(logP),1),repmat(mw,numel(logP),1), ...
    logP,logPrediction,psbbb,'VariableNames', ...
    {'Model','MW_g_per_mol','logP','Predicted_log10_parameter','PS_BBB_L_per_h'});
end

function [handles,labels] = modelLegend(ax,mwLines,lineStyles,colors,lineWidth,grumettoWidth)
names = ["Smith","Luco","Tsinman"];
handles = gobjects(10,1); labels = strings(10,1); index = 1;
for m = 1:numel(names)
    for i = 1:numel(mwLines)
        handles(index) = plot(ax,nan,nan,'Color',colors.(names(m)), ...
            'LineStyle',lineStyles{i},'LineWidth',lineWidth);
        labels(index) = names(m)+" MW "+mwLines(i); index = index+1;
    end
end
handles(index) = plot(ax,nan,nan,'Color',colors.Grumetto, ...
    'LineStyle','-','LineWidth',grumettoWidth);
labels(index) = "Grumetto";
labels = cellstr(labels);
end
