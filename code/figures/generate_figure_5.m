%% Generate Figure 5: GMR and GMFE comparison

clear; clc; close all;
rng(12,'twister')

%% ===============================================================
% GLOBAL STYLE
% ===============================================================

FIG_WIDTH  = 17.8;   % cm
FIG_HEIGHT = 11;    % cm

AXIS_FONT = 8;
Y_LABEL_OFFSET = 3;
Y_LABEL_SIZE = 8;
AXIS_LINE_WIDTH = 1.2;

LEGEND_FONT = 8;

MARKER_SIZE = 5;
MARKER_SIZE_BIG = 10;

LINE_WIDTH_MARKER = 1.5;
LINE_WIDTH_MARKER_BIG = 1.7;

TITLE_FONT = 10;

BAR_LINE_WIDTH = 2;
SD_LINE_WIDTH = 1.5;

%% ===============================================================
% Load data
% ===============================================================

scriptDir = fileparts(mfilename('fullpath'));
repoRoot = fileparts(fileparts(scriptDir));
dataFile = fullfile(repoRoot, 'results', 'derived', 'main_figures', ...
    'figure_5', 'fold_error_summary.csv');
outputDir = fullfile(repoRoot, 'results', 'figures', 'main');

if ~isfile(dataFile)
    error('Required fold-error summary was not found: %s', dataFile);
end
if ~exist(outputDir, 'dir')
    mkdir(outputDir);
end

T = readtable(dataFile);
T.Compound(strcmp(T.Compound,'Indometacin')) = {'Indomethacin'};

isOpt = strcmp(T.Source,'optimized');
T.isOpt = isOpt;
T = sortrows(T, {'Compound','isOpt'}, {'ascend','ascend'});
T.isOpt = [];

compounds = {'Ethanol','Diazepam','Paracetamol','Ibuprofen','Indomethacin','Mannitol'};
sources   = {'optimized','perfusion-limited','in situ perfusion',...
             'in vitro cell-based','in vitro PAMPA',...
             'in silico perfusion','in silico PAMPA'};

xpos_comp = containers.Map(compounds,1:numel(compounds));
xpos_src  = containers.Map(sources,1:numel(sources));


%% ===============================================================
% Colors 
% ===============================================================

faceMap = containers.Map(...
    {'in situ perfusion','in vitro cell-based','in vitro PAMPA', ...
     'in silico perfusion','in silico PAMPA','perfusion-limited'}, ...
    {[0.9255 0.0039 0.0039],[0 0.4235 0.9451],[0 0.4745 0],'none','none','none'});

edgeMap = containers.Map(...
    {'in situ perfusion','in vitro cell-based','in vitro PAMPA', ...
     'in silico perfusion','in silico PAMPA','perfusion-limited'}, ...
    {[0.9255 0.0039 0.0039], [0 0.4235 0.9451], [0 0.4745 0],[0.9255 0.0039 0.0039],[0 0.4745 0],'k'});

%% ===============================================================
% Fold transform
% ===============================================================

compress = 0.1;

fold_transform = @(y) arrayfun(@(v) local_transform(v,compress), y);

%% ===============================================================
% Prepare data
% ===============================================================

N = height(T);

X_comp = zeros(N,1);
for i=1:N
    X_comp(i) = xpos_comp(T.Compound{i});
end

GMR_raw = T.GMR;
GMR_transformed = fold_transform(GMR_raw);

%% ===============================================================
% CREATE FIGURE
% ===============================================================

figure
set(gcf,'Units','centimeters','Position',[2 2 FIG_WIDTH FIG_HEIGHT])

t = tiledlayout(1,2,'Padding','compact','TileSpacing','compact');

%% ===============================================================
% LEFT PANEL
% ===============================================================
axLeft = nexttile;
hold(axLeft,'on'); box(axLeft,'on');

xfill = [0.5 numel(compounds)+0.5 numel(compounds)+0.5 0.5];

fill(xfill,fold_transform([1/5 1/5 1/3 1/3]),[1 0 0],'FaceAlpha',0.15,'EdgeColor','none');
fill(xfill,fold_transform([1/3 1/3 0.5 0.5]),[1 0.8 0],'FaceAlpha',0.15,'EdgeColor','none');
fill(xfill,fold_transform([0.5 0.5 1 1]),[0 0.8 0],'FaceAlpha',0.15,'EdgeColor','none');
fill(xfill,fold_transform([1 1 2 2]),[0 0.8 0],'FaceAlpha',0.15,'EdgeColor','none');
fill(xfill,fold_transform([2 2 3 3]),[1 0.8 0],'FaceAlpha',0.15,'EdgeColor','none');
fill(xfill,fold_transform([3 3 5 5]),[1 0 0],'FaceAlpha',0.15,'EdgeColor','none');

yline(fold_transform(1),'k--','LineWidth',1.2);

for i=1:N
    src = T.Source{i};
    rnd = rand;
    Xi = X_comp(i) + (rnd*0.2-0.1);
    Yi = GMR_transformed(i);

    if strcmp(src,'optimized')
        Xi = Xi - (rnd*0.2-0.1);
        plot(Xi,Yi,'kx','MarkerSize',MARKER_SIZE_BIG,'LineWidth',LINE_WIDTH_MARKER_BIG);

    elseif strcmp(src,'perfusion-limited')
        Xi = Xi - (rnd*0.2-0.1);
        plot(Xi,Yi,'ks','MarkerSize',MARKER_SIZE_BIG,...
             'MarkerFaceColor','none','LineWidth',LINE_WIDTH_MARKER_BIG);

    elseif isKey(faceMap,src)
        faceCol = faceMap(src);
        edgeCol = edgeMap(src);

        if strcmp(faceCol,'none')
            plot(Xi,Yi,'o','MarkerFaceColor','none',...
                'MarkerEdgeColor',edgeCol,...
                'MarkerSize',MARKER_SIZE,'LineWidth',LINE_WIDTH_MARKER);
        else
            plot(Xi,Yi,'o','MarkerFaceColor',faceCol,...
                'MarkerEdgeColor',edgeCol,...
                'MarkerSize',MARKER_SIZE, 'LineWidth',LINE_WIDTH_MARKER);
        end
    end
end

xlim([0.5 numel(compounds)+0.5]);
ylim(fold_transform([1/100 100]));
set(gca,'XTick',1:numel(compounds),'XTickLabel',compounds,'FontSize',AXIS_FONT, 'LineWidth', AXIS_LINE_WIDTH);

tickVals = [1/50 1/5 1/3 0.5 1 2 3 5 50];
yticks(fold_transform(tickVals));

yticklabels({'1/50','1/5','1/3','1/2',...
             '1','2','3','5','50'});

ax = gca;
ax.YAxis.TickLabelGapOffset = Y_LABEL_OFFSET;

ylabel('GMR', 'FontSize', Y_LABEL_SIZE);

hPanelA = text(axLeft,0,0,'A', ...
    'HorizontalAlignment','left', ...
    'VerticalAlignment','bottom', ...
    'FontSize',TITLE_FONT + 4, ...
    'FontWeight','bold', ...
    'Clipping','off');
hPanelA.Units = 'normalized';
hPanelA.Position = [-0.05 1.003 0];

%% ===============================================================
% LEGEND 
% ===============================================================

sourceHandles = gobjects(numel(sources),1);

for s = 1:numel(sources)
    src = sources{s};

    if strcmp(src,'optimized')
        sourceHandles(s) = plot(nan,nan,'kx','MarkerSize',MARKER_SIZE_BIG,'LineWidth',LINE_WIDTH_MARKER_BIG);

    elseif strcmp(src,'perfusion-limited')
        sourceHandles(s) = plot(nan,nan,'ks','MarkerFaceColor','none','MarkerSize',MARKER_SIZE_BIG,'LineWidth',LINE_WIDTH_MARKER_BIG);

    elseif isKey(faceMap,src)
        faceCol = faceMap(src);
        edgeCol = edgeMap(src);

        if strcmp(faceCol,'none')
            sourceHandles(s) = plot(nan,nan,'o','MarkerFaceColor','none',...
                'MarkerEdgeColor',edgeCol,'MarkerSize',MARKER_SIZE,'LineWidth',LINE_WIDTH_MARKER);
        else
            sourceHandles(s) = plot(nan,nan,'o','MarkerFaceColor',faceCol,...
                'MarkerEdgeColor',edgeCol,'MarkerSize',MARKER_SIZE,'LineWidth',LINE_WIDTH_MARKER);
        end
    end
end

lgd = legend(sourceHandles, sources,...
    'Orientation','horizontal',...
    'Location','northoutside',...
    'Box','on',...
    'FontSize',LEGEND_FONT);

lgd.Layout.Tile = 'north';
lgd.NumColumns = 4;

%% ===============================================================
% RIGHT PANEL 
% ===============================================================

axRight = nexttile;
hold(axRight,'on'); box(axRight,'on');

GMFE_raw = T.GMFE;

mean_GMFE = zeros(numel(sources),1);
sd_GMFE   = zeros(numel(sources),1);

for s = 1:numel(sources)
    idx = strcmp(T.Source, sources{s});
    vals = GMFE_raw(idx);

    mean_GMFE(s) = mean(vals,'omitnan');
    sd_GMFE(s)   = std(vals,'omitnan');
end

[mean_GMFE_sorted, sortIdx] = sort(mean_GMFE,'ascend');
sd_GMFE_sorted = sd_GMFE(sortIdx);
sources_sorted = sources(sortIdx);

Y_mean = fold_transform(mean_GMFE_sorted);
Y_sd_upper = fold_transform(mean_GMFE_sorted + sd_GMFE_sorted);
Y_sd_lower = fold_transform(max(mean_GMFE_sorted - sd_GMFE_sorted,1));

xfill = [0.5 numel(sources_sorted)+0.5 numel(sources_sorted)+0.5 0.5];

fill(xfill,fold_transform([1 1 2 2]),[0 0.8 0],'FaceAlpha',0.15,'EdgeColor','none');
fill(xfill,fold_transform([2 2 3 3]),[1 0.8 0],'FaceAlpha',0.15,'EdgeColor','none');
fill(xfill,fold_transform([3 3 5 5]),[1 0 0],'FaceAlpha',0.15,'EdgeColor','none');

for i = 1:numel(sources_sorted)
    src = sources_sorted{i};
    y   = Y_mean(i);

    if strcmp(src,'optimized')
        bar(i,y,'FaceColor','k','EdgeColor','k','LineWidth',BAR_LINE_WIDTH);

    elseif strcmp(src,'perfusion-limited')
        bar(i,y,'FaceColor','none','EdgeColor','k','LineWidth',BAR_LINE_WIDTH);

    elseif isKey(faceMap,src)
        faceCol = faceMap(src);
        edgeCol = edgeMap(src);

        if strcmp(faceCol,'none')
            bar(i,y,'FaceColor','none','EdgeColor',edgeCol,'LineWidth',BAR_LINE_WIDTH);
        else
            bar(i,y,'FaceColor',faceCol,'EdgeColor',edgeCol,'LineWidth',BAR_LINE_WIDTH);
        end
    end
end

for i = 1:numel(sources_sorted)
    plot([i i],[Y_sd_lower(i) Y_sd_upper(i)],'k','LineWidth',SD_LINE_WIDTH);
end

xlim([0.5 numel(sources_sorted)+0.5]);
ylim(fold_transform([1 100]));

set(gca,'XTick',1:numel(sources_sorted),...
    'XTickLabel',sources_sorted,...
    'XTickLabelRotation',45,...
    'FontSize',AXIS_FONT, ...
    'LineWidth', AXIS_LINE_WIDTH);

tickVals = [1 2 3 5 50];
yticks(fold_transform(tickVals));

yticklabels({'1','2','3','5','50'});

ax = gca;
ax.YAxis.TickLabelGapOffset = Y_LABEL_OFFSET;

ylabel('GMFE', 'FontSize',  Y_LABEL_SIZE);

hPanelB = text(axRight,0,0,'B', ...
    'HorizontalAlignment','left', ...
    'VerticalAlignment','bottom', ...
    'FontSize',TITLE_FONT + 4, ...
    'FontWeight','bold', ...
    'Clipping','off');
hPanelB.Units = 'normalized';
hPanelB.Position = [-0.05 1.003 0];

%% ===============================================================
% EXPORT PDF
% ===============================================================

set(gcf,'PaperUnits','centimeters')
set(gcf,'PaperPosition',[0 0 FIG_WIDTH FIG_HEIGHT])
set(gcf,'PaperSize',[FIG_WIDTH FIG_HEIGHT])

pdfFile = fullfile(outputDir, 'Figure_5.pdf');
pngFile = fullfile(outputDir, 'Figure_5.png');
print(gcf, pdfFile, '-dpdf', '-vector')
exportgraphics(gcf, pngFile, 'Resolution', 600)

fprintf('Saved %s\n', pdfFile);
fprintf('Saved %s\n', pngFile);

%% ===============================================================
% Local functions
% ===============================================================

function yt = local_transform(v,c)
    L  = log10(v);
    L5 = log10(5);
    if abs(L) <= L5
        yt = L;
    else
        yt = sign(L) * ( L5 + c*(abs(L)-L5) );
    end
end
