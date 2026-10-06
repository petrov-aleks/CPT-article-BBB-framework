%% Generate Figure S15: predictive performance of alternative models

clear; clc; close all; rng(12,'twister');
scriptDir=fileparts(mfilename('fullpath'));
repoRoot=fileparts(fileparts(scriptDir));
inputDir=fullfile(repoRoot,'results','derived','supplementary_figures','figure_s15');
compoundFile=fullfile(inputDir,'figure_s15_compound_summary.csv');
modelFile=fullfile(inputDir,'figure_s15_model_summary.csv');
outputDir=fullfile(repoRoot,'results','figures','supplementary');
assert(isfile(compoundFile) && isfile(modelFile), ...
    'Run code/analysis/calculate_figure_s15_summary.m first.');
if ~exist(outputDir,'dir'), mkdir(outputDir); end
compoundData=readtable(compoundFile,'TextType','string');
modelData=readtable(modelFile,'TextType','string');

compounds=["Ethanol","Diazepam","Paracetamol", ...
    "Ibuprofen","Indomethacin","Mannitol"];
perfusionSources=["in silico perfusion (Smith et al. 2024)", ...
    "in silico perfusion (Luco et al. 2006)"];
pampaSources=["in silico PAMPA (Tsinman et al. 2011)", ...
    "in silico PAMPA (Grumetto et al. 2016)"];
colors=containers.Map(cellstr([perfusionSources pampaSources]), ...
    {[0.9255 0.0039 0.0039],[0.6510 0.0039 0.4745], ...
     [0 0.4745 0],[0.9961 0.7059 0.3843]});

fig=figure('Color','w','Units','centimeters','Position',[2 2 17.8 18]);
layout=tiledlayout(fig,2,2,'Padding','compact','TileSpacing','compact');
plotGMR(nexttile(layout,1),compoundData,compounds,perfusionSources,colors);
plotGMFE(nexttile(layout,2),modelData,perfusionSources,colors);
plotGMR(nexttile(layout,3),compoundData,compounds,pampaSources,colors);
plotGMFE(nexttile(layout,4),modelData,pampaSources,colors);

drawnow;
allAxes=findall(fig,'Type','axes');
for i=1:numel(allAxes)
    allAxes(i).Toolbar.Visible='off';
end

pdfFile=fullfile(outputDir,'Figure_S15.pdf');
pngFile=fullfile(outputDir,'Figure_S15.png');
exportgraphics(fig,pdfFile,'ContentType','vector');
exportgraphics(fig,pngFile,'Resolution',600);
fprintf('Saved %s\nSaved %s\n',pdfFile,pngFile);

function plotGMR(ax,data,compounds,sources,colors)
hold(ax,'on'); box(ax,'on'); addGMRZones(ax,numel(compounds));
yline(ax,foldTransform(1),'k--','LineWidth',1.2);
% Preserve the row traversal used by the original plotting workflow so that
% rng(12) assigns exactly the same horizontal jitter to every point.
legacyCompoundOrder=["Ethanol","Paracetamol","Mannitol", ...
    "Ibuprofen","Diazepam","Indomethacin"];
data=data(ismember(data.Source,sources),:);
data.LegacyCompoundRank=zeros(height(data),1);
for i=1:height(data)
 data.LegacyCompoundRank(i)=find(legacyCompoundOrder==data.Compound(i),1);
end
data=sortrows(data,{'LegacyCompoundRank','PS_BBB_L_per_h'},{'ascend','ascend'});
for i=1:height(data)
 x=find(compounds==data.Compound(i),1)+(rand*0.2-0.1);
 plot(ax,x,foldTransform(data.GMR(i)),'o','MarkerFaceColor','none', ...
     'MarkerEdgeColor',colors(char(data.Source(i))),'MarkerSize',6,'LineWidth',1.5);
end
xlim(ax,[0.5 numel(compounds)+0.5]); ylim(ax,foldTransform([1/100 100]));
set(ax,'XTick',1:numel(compounds),'XTickLabel',cellstr(compounds), ...
    'FontSize',8,'LineWidth',1.2);
ticks=[1/50 1/5 1/3 0.5 1 2 3 5 50]; yticks(ax,foldTransform(ticks));
yticklabels(ax,{'1/50','1/5','1/3','1/2','1','2','3','5','50'});
ylabel(ax,'GMR');
handles=gobjects(numel(sources),1);
for i=1:numel(sources)
 handles(i)=plot(ax,nan,nan,'o','MarkerFaceColor','none', ...
     'MarkerEdgeColor',colors(char(sources(i))),'MarkerSize',6,'LineWidth',1.7);
end
lgd=legend(ax,handles,cellstr(sources),'Orientation','horizontal', ...
    'Location','northoutside','Box','on','FontSize',6);
lgd.Layout.Tile='north';
end

function plotGMFE(ax,data,sources,colors)
hold(ax,'on'); box(ax,'on'); addGMFEZones(ax,numel(sources));
for i=1:numel(sources)
 row=data(data.Source==sources(i),:);
 assert(height(row)==1,'Missing model summary for %s.',sources(i));
 meanValue=row.Mean_GMFE; sdValue=row.SD_GMFE;
 bar(ax,i,foldTransform(meanValue),'FaceColor','none', ...
     'EdgeColor',colors(char(sources(i))),'LineWidth',2);
 lower=max(meanValue-sdValue,1); upper=meanValue+sdValue;
 plot(ax,[i i],[foldTransform(lower) foldTransform(upper)],'k','LineWidth',1.5);
end
xlim(ax,[0.5 numel(sources)+0.5]); ylim(ax,foldTransform([1 100]));
labels=cell(numel(sources),1);
for i=1:numel(sources)
    parts=regexp(char(sources(i)),'^(.*) \((.*)\)$','tokens','once');
    labels{i}=sprintf('%s\\newline(%s)',parts{1},parts{2});
end
set(ax,'XTick',1:numel(sources),'XTickLabel',labels, ...
    'XTickLabelRotation',0,'TickLabelInterpreter','tex','FontSize',8,'LineWidth',1.2);
ticks=[1 2 3 5 50]; yticks(ax,foldTransform(ticks));
yticklabels(ax,{'1','2','3','5','50'}); ylabel(ax,'GMFE');
end

function addGMRZones(ax,n)
x=[0.5 n+0.5 n+0.5 0.5];
addZone(ax,x,1/5,1/3,[1 0 0]); addZone(ax,x,1/3,0.5,[1 0.8 0]);
addZone(ax,x,0.5,1,[0 0.8 0]); addZone(ax,x,1,2,[0 0.8 0]);
addZone(ax,x,2,3,[1 0.8 0]); addZone(ax,x,3,5,[1 0 0]);
end

function addGMFEZones(ax,n)
x=[0.5 n+0.5 n+0.5 0.5];
addZone(ax,x,1,2,[0 0.8 0]); addZone(ax,x,2,3,[1 0.8 0]); addZone(ax,x,3,5,[1 0 0]);
end

function addZone(ax,x,low,high,color)
fill(ax,x,foldTransform([low low high high]),color, ...
    'FaceAlpha',0.15,'EdgeColor','none','HandleVisibility','off');
end

function transformed=foldTransform(values)
compression=0.1; boundary=log10(5); logs=log10(values);
transformed=logs;
outside=abs(logs)>boundary;
transformed(outside)=sign(logs(outside)).*(boundary+compression.*(abs(logs(outside))-boundary));
end
