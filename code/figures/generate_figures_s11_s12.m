%% Generate Figures S11 and S12 from the fold-change summary data

clear; clc; rng(12,'twister');
scriptDir=fileparts(mfilename('fullpath'));
repoRoot=fileparts(fileparts(scriptDir));
inputFile=fullfile(repoRoot,'results','summary_data','figures_s11_s12_fold_changes.csv');
if ~isfile(inputFile)
    error(['Summary data were not found. Run the two scripts in code/analysis ' ...
        'for Figures S11 and S12 first.']);
end
data=readtable(inputFile,'TextType','string');
outputDir=fullfile(repoRoot,'results','figures','supplementary');
if ~exist(outputDir,'dir'), mkdir(outputDir); end

makeFigure(data,"brain ECF","Figure_S11",outputDir);
makeFigure(data,"scsf","Figure_S12",outputDir);

function makeFigure(data,site,outputName,outputDir)
FIG_WIDTH=17.8; FIG_HEIGHT=18; AXIS_FONT=8; AXIS_LINE_WIDTH=1.2;
MARKER_SIZE=6; MARKER_SIZE_BIG=12; MARKER_WIDTH=1.5; BIG_WIDTH=1.7;
BAR_WIDTH=2; SD_WIDTH=1.5; TITLE_FONT=12; LEGEND_FONT=7;
compounds=["Ethanol","Diazepam","Paracetamol","Ibuprofen","Indomethacin","Mannitol"];
sources=["optimized","perfusion-limited","in situ perfusion", ...
    "in vitro cell-based","in vitro PAMPA","in silico perfusion","in silico PAMPA"];
metrics=["Cmax","AUC"];
T=data(strcmpi(data.Site,site),:);

fig=figure('Units','centimeters','Position',[2 2 FIG_WIDTH FIG_HEIGHT]);
tl=tiledlayout(fig,2,2,'Padding','compact','TileSpacing','compact');
for m=1:2
    D=T(T.Metric==metrics(m),:);
    ax=nexttile(tl); hold(ax,'on'); box(ax,'on');
    addRatioBands(ax,numel(compounds)); yline(ax,foldTransform(1),'k--','LineWidth',1.2);
    for i=1:height(D)
        x=find(compounds==D.Compound(i));
        if isempty(x), continue; end
        if any(D.Source(i)==["optimized","perfusion-limited"]), jitter=0;
        else, jitter=0.2*(rand-0.5); end
        drawMarker(ax,x+jitter,foldTransform(D.Ratio(i)),D.Source(i), ...
            MARKER_SIZE,MARKER_SIZE_BIG,MARKER_WIDTH,BIG_WIDTH);
    end
    xlim(ax,[0.5 numel(compounds)+0.5]); ylim(ax,foldTransform([1/100 100]));
    set(ax,'XTick',1:numel(compounds),'XTickLabel',compounds,'FontSize',AXIS_FONT, ...
        'LineWidth',AXIS_LINE_WIDTH);
    yticks(ax,foldTransform([1/50 1/5 1/3 0.5 1 2 3 5 50]));
    yticklabels(ax,{'1/50','1/5','1/3','1/2','1','2','3','5','50'});
    ylabel(ax,'R'); title(ax,metricTitle(metrics(m)),'Interpreter','tex','FontSize',TITLE_FONT);

    ax=nexttile(tl); hold(ax,'on'); box(ax,'on');
    means=nan(numel(sources),1); deviations=nan(numel(sources),1);
    for s=1:numel(sources)
        values=D.FoldError(D.Source==sources(s));
        means(s)=mean(values,'omitnan'); deviations(s)=std(values,'omitnan');
    end
    [means,order]=sort(means,'ascend','MissingPlacement','last');
    deviations=deviations(order); sortedSources=sources(order);
    addFoldErrorBands(ax,numel(sources));
    for s=1:numel(sortedSources)
        if isnan(means(s)), continue; end
        drawBar(ax,s,foldTransform(means(s)),sortedSources(s),BAR_WIDTH);
        lo=foldTransform(max(means(s)-deviations(s),1));
        hi=foldTransform(means(s)+deviations(s));
        plot(ax,[s s],[lo hi],'k','LineWidth',SD_WIDTH);
    end
    xlim(ax,[0.5 numel(sources)+0.5]); ylim(ax,foldTransform([1 100]));
    set(ax,'XTick',1:numel(sources),'XTickLabel',sortedSources,'XTickLabelRotation',45, ...
        'FontSize',AXIS_FONT,'LineWidth',AXIS_LINE_WIDTH);
    yticks(ax,foldTransform([1 2 3 5 50])); yticklabels(ax,{'1','2','3','5','50'});
    ylabel(ax,'FE'); title(ax,metricTitle(metrics(m)),'Interpreter','tex','FontSize',TITLE_FONT);
end

handles=gobjects(numel(sources),1);
for s=1:numel(sources)
    handles(s)=legendMarker(sources(s),MARKER_SIZE,MARKER_SIZE_BIG,MARKER_WIDTH,BIG_WIDTH);
end
lgd=legend(handles,sources,'Orientation','horizontal','Location','northoutside', ...
    'Box','on','FontSize',LEGEND_FONT); lgd.Layout.Tile='north';
exportgraphics(fig,fullfile(outputDir,outputName+'.pdf'),'ContentType','vector');
exportgraphics(fig,fullfile(outputDir,outputName+'.png'),'Resolution',600);
fprintf('%s saved in %s\n',outputName,outputDir);
close(fig);
end

function addRatioBands(ax,n)
x=[0.5 n+0.5 n+0.5 0.5];
fill(ax,x,foldTransform([1/5 1/5 1/3 1/3]),[1 0 0],'FaceAlpha',.15,'EdgeColor','none');
fill(ax,x,foldTransform([1/3 1/3 .5 .5]),[1 .8 0],'FaceAlpha',.15,'EdgeColor','none');
fill(ax,x,foldTransform([.5 .5 1 1]),[0 .8 0],'FaceAlpha',.15,'EdgeColor','none');
fill(ax,x,foldTransform([1 1 2 2]),[0 .8 0],'FaceAlpha',.15,'EdgeColor','none');
fill(ax,x,foldTransform([2 2 3 3]),[1 .8 0],'FaceAlpha',.15,'EdgeColor','none');
fill(ax,x,foldTransform([3 3 5 5]),[1 0 0],'FaceAlpha',.15,'EdgeColor','none');
end

function addFoldErrorBands(ax,n)
x=[0.5 n+0.5 n+0.5 0.5];
fill(ax,x,foldTransform([1 1 2 2]),[0 .8 0],'FaceAlpha',.15,'EdgeColor','none');
fill(ax,x,foldTransform([2 2 3 3]),[1 .8 0],'FaceAlpha',.15,'EdgeColor','none');
fill(ax,x,foldTransform([3 3 5 5]),[1 0 0],'FaceAlpha',.15,'EdgeColor','none');
end

function drawMarker(ax,x,y,source,ms,msBig,lw,lwBig)
[face,edge]=sourceColors(source);
if source=="optimized"
    plot(ax,x,y,'kx','MarkerSize',msBig,'LineWidth',lwBig);
elseif source=="perfusion-limited"
    plot(ax,x,y,'ks','MarkerFaceColor','none','MarkerSize',msBig,'LineWidth',lwBig);
else
    plot(ax,x,y,'o','MarkerFaceColor',face,'MarkerEdgeColor',edge,'MarkerSize',ms,'LineWidth',lw);
end
end

function drawBar(ax,x,y,source,lw)
[face,edge]=sourceColors(source);
if source=="optimized"
    face='k'; edge='k';
elseif source=="perfusion-limited"
    face='none'; edge='k';
end
bar(ax,x,y,'FaceColor',face,'EdgeColor',edge,'LineWidth',lw);
end

function h=legendMarker(source,ms,msBig,lw,lwBig)
[face,edge]=sourceColors(source);
if source=="optimized"
    h=plot(nan,nan,'kx','MarkerSize',msBig,'LineWidth',lwBig);
elseif source=="perfusion-limited"
    h=plot(nan,nan,'ks','MarkerFaceColor','none','MarkerSize',msBig,'LineWidth',lwBig);
else
    h=plot(nan,nan,'o','MarkerFaceColor',face,'MarkerEdgeColor',edge,'MarkerSize',ms,'LineWidth',lw);
end
end

function [face,edge]=sourceColors(source)
red=[.9255 .0039 .0039]; blue=[0 .4235 .9451]; green=[0 .4745 0];
switch source
 case "in situ perfusion", face=red; edge=red;
 case "in vitro cell-based", face=blue; edge=blue;
 case "in vitro PAMPA", face=green; edge=green;
 case "in silico perfusion", face='none'; edge=red;
 case "in silico PAMPA", face='none'; edge=green;
 otherwise, face='none'; edge='k';
end
end

function y=foldTransform(x)
c=.1; L=log10(x); L5=log10(5); y=L;
idx=abs(L)>L5; y(idx)=sign(L(idx)).*(L5+c*(abs(L(idx))-L5));
end

function value=metricTitle(metric)
if metric=="Cmax", value='C_{max}'; else, value='AUC'; end
end
