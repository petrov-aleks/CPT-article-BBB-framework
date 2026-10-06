%% Generate Figure 3: intrinsic and effective BBB permeability

clear; clc; close all;

% -------------------------------------------------------------------------
% GLOBAL FIGURE STYLE PARAMETERS (CPT:PSP journal compliant)
% -------------------------------------------------------------------------

FIG_WIDTH  = 17.8;   % cm (double column)
FIG_HEIGHT = 16.5;    % cm

AXIS_LINE_WIDTH = 1.2;
BOX_LINE_WIDTH  = 1.5;

X_TICK_FONT = 10;
Y_TICK_FONT = 10;

AXIS_LABEL_FONT = 10;
TITLE_FONT = 12;
LEGEND_FONT = 10;

FU_FN_FONT = 9;
ZONE_FONT = 10;
FORMULA_FONT = 10;

POINT_SIZE = 45;
POINT_LINE_WIDTH = 1.5;

% -------------------------------------------------------------------------
% Load data
% -------------------------------------------------------------------------

scriptDir = fileparts(mfilename('fullpath'));
repoRoot = fileparts(fileparts(scriptDir));
tableFile = fullfile(repoRoot, 'results', 'tables', ...
    'figure_3_permeability.csv');
outputDir = fullfile(repoRoot, 'results', 'figures', 'main');

if ~isfile(tableFile)
    error(['Required input was not found: %s\n' ...
        'Run code/tables/generate_table_1.m first.'], tableFile);
end
if ~exist(outputDir, 'dir')
    mkdir(outputDir);
end

T = readtable(tableFile, 'TextType', 'string');
T.Compound = string(T.Compound);

% -------------------------------------------------------------------------
% Define compounds
% -------------------------------------------------------------------------

compounds = ["Ethanol","Paracetamol","Mannitol", ...
             "Ibuprofen","Diazepam","Indomethacin"];

% -------------------------------------------------------------------------
% Scaling parameters
% -------------------------------------------------------------------------

params = table( ...
    ["Ethanol";"Paracetamol";"Mannitol"; ...
     "Ibuprofen";"Diazepam";"Indomethacin"], ...
    [1;0.86;1;0.01;0.013;0.01], ...
    [1;1;1;0.0079;1;0.0013], ...
    'VariableNames',{'Compound','fuP','fn'});

% -------------------------------------------------------------------------
% Extract data
% -------------------------------------------------------------------------

PS_raw_all    = [];
PS_scaled_all = [];
group_labels  = strings(0,1);
Sources_all   = strings(0,1);

median_insitu_raw    = zeros(numel(compounds),1);
median_insitu_scaled = zeros(numel(compounds),1);

for i = 1:numel(compounds)

    compound = compounds(i);

    idx_data = strcmpi(T.Compound,compound);
    PS_raw = T.PS_BBB_L_per_h(idx_data);
    Sources = T.Source(idx_data);

    idx_param = strcmpi(params.Compound,compound);
    fuP = params.fuP(idx_param);
    fn  = params.fn(idx_param);

    PS_scaled = fuP * fn .* PS_raw;

    median_insitu_raw(i)    = median(PS_raw(strcmpi(Sources,'in situ perfusion')),'omitnan');
    median_insitu_scaled(i) = median(PS_scaled(strcmpi(Sources,'in situ perfusion')),'omitnan');

    PS_raw_all    = [PS_raw_all; PS_raw];
    PS_scaled_all = [PS_scaled_all; PS_scaled];
    group_labels  = [group_labels; repmat(compound,numel(PS_raw),1)];
    Sources_all   = [Sources_all; Sources];

end

[~, idx_raw_order]    = sort(median_insitu_raw,'descend');
[~, idx_scaled_order] = sort(median_insitu_scaled,'descend');

order_raw    = compounds(idx_raw_order);
order_scaled = compounds(idx_scaled_order);

group_cat_raw    = categorical(group_labels, order_raw, 'Ordinal',true);
group_cat_scaled = categorical(group_labels, order_scaled, 'Ordinal',true);

[~, cat_map] = ismember(order_raw, order_scaled);

% -------------------------------------------------------------------------
% Create figure
% -------------------------------------------------------------------------

figure
set(gcf,'Units','centimeters','Position',[2 2 FIG_WIDTH FIG_HEIGHT])

tiledlayout(2,1,'Padding','compact','TileSpacing','compact');

rng(12,'twister')

% -------------------------------------------------------------------------
% PANEL A
% -------------------------------------------------------------------------

nexttile

h1 = boxplot(PS_raw_all, group_cat_raw,'Colors','k','Symbol','');
set(h1,'LineWidth',BOX_LINE_WIDTH);

set(gca,'YScale','log','LineWidth',AXIS_LINE_WIDTH)

ax = gca;
ax.XAxis.FontSize = X_TICK_FONT;
ax.YAxis.FontSize = Y_TICK_FONT;

ax.YMinorTick = 'off';

ylim([1e-6 1e6])
yticks(10.^([-6 -4 -2 0 2 4 6]))

ylabel('PS_{BBB} (L/h)','FontSize',AXIS_LABEL_FONT)

hold on

x_raw = double(group_cat_raw);
jitter = (rand(numel(PS_raw_all),1)-0.5)*0.25;
xj_raw = x_raw + jitter;

for i = 1:numel(PS_raw_all)

[c, filled] = mapSourceToStyle(Sources_all(i));

if filled
scatter(xj_raw(i),PS_raw_all(i),POINT_SIZE,c,'filled')
else
scatter(xj_raw(i),PS_raw_all(i),POINT_SIZE,c,'LineWidth',POINT_LINE_WIDTH)
end

end

title('A. Intrinsic BBB Permeability-Surface Area Product','FontSize',TITLE_FONT,'FontWeight','bold')

% fu fn

xTicks = 1:numel(order_raw);
yLimits = ylim;

yPosText = 10^( log10(yLimits(1)) + 0.1*(log10(yLimits(2))-log10(yLimits(1))) );

for i = 1:numel(order_raw)

compound = order_raw(i);
idx_param = strcmpi(params.Compound,compound);

fu_val = params.fuP(idx_param);
fn_val = params.fn(idx_param);

% text(xTicks(i),yPosText,...
% sprintf('f_u = %.1f%%\nf_n = %.1f%%',fu_val*100,fn_val*100),...
% 'HorizontalAlignment','center',...
% 'FontSize',FU_FN_FONT,'FontWeight','bold');

text(xTicks(i),yPosText,...
    sprintf('fu_{brp} = %.1f%%\nfn_{brp} = %.1f%%',fu_val*100,fn_val*100),...
    'HorizontalAlignment','center',...
    'FontSize',FU_FN_FONT,'FontWeight','normal');


end

% legend

hL1 = scatter(NaN,NaN,POINT_SIZE,[0.9255 0.0039 0.0039],'filled');
hL2 = scatter(NaN,NaN,POINT_SIZE,[0 0.4235 0.9451],'filled');
hL3 = scatter(NaN,NaN,POINT_SIZE,[0 0.4745 0],'filled');
hL4 = scatter(NaN,NaN,POINT_SIZE,[0.9255 0.0039 0.0039],'o','LineWidth',POINT_LINE_WIDTH);
hL5 = scatter(NaN,NaN,POINT_SIZE,[0 0.4745 0],'o','LineWidth',POINT_LINE_WIDTH);

legend([hL1 hL2 hL3 hL4 hL5],...
{'in situ perfusion','in vitro cell-based','in vitro PAMPA',...
'in silico perfusion','in silico PAMPA'},...
'Location','northoutside','Orientation','horizontal','Box','on','FontSize',LEGEND_FONT)

% -------------------------------------------------------------------------
% PANEL B
% -------------------------------------------------------------------------

nexttile

h2 = boxplot(PS_scaled_all, group_cat_scaled,'Colors','k','Symbol','');
set(h2,'LineWidth',BOX_LINE_WIDTH);

set(gca,'YScale','log','LineWidth',AXIS_LINE_WIDTH)

ax = gca;
ax.XAxis.FontSize = X_TICK_FONT;
ax.YAxis.FontSize = Y_TICK_FONT;

ax.YMinorTick = 'off';

ylim([1e-6 1e6])
yticks(10.^([-6 -4 -2 0 2 4 6]))

ylabel('fu_{brp} × fn_{brp} × PS_{BBB} (L/h)','FontSize',AXIS_LABEL_FONT)


hold on

x_scaled = double(group_cat_scaled);
jitter_scaled = jitter;

for k = 1:numel(order_raw)

    idx = strcmp(group_labels, order_raw(k));
    jitter_scaled(idx) = jitter(strcmp(group_labels, order_scaled(cat_map(k))));

end

xj_scaled = x_scaled + jitter_scaled;

for i = 1:numel(PS_scaled_all)

[c, filled] = mapSourceToStyle(Sources_all(i));

if filled
scatter(xj_scaled(i),PS_scaled_all(i),POINT_SIZE,c,'filled')
else
scatter(xj_scaled(i),PS_scaled_all(i),POINT_SIZE,c,'LineWidth',POINT_LINE_WIDTH)
end

end

title('B. Effective BBB Permeability-Surface Area Product','FontSize',TITLE_FONT,'FontWeight','bold')

% -------------------------------------------------------------------------
% Zones
% -------------------------------------------------------------------------

PsPermLimit = 5;
PsPerfLimit = 140;

permColor = [1 0.5 0];
perfColor = [0.2 0.2 0.2];

xLimits = xlim;
yLimits = ylim;

patch([xLimits(1) xLimits(2) xLimits(2) xLimits(1)],...
[yLimits(1) yLimits(1) PsPermLimit PsPermLimit],...
permColor,'FaceAlpha',0.08,'EdgeColor','none')

patch([xLimits(1) xLimits(2) xLimits(2) xLimits(1)],...
[PsPerfLimit PsPerfLimit yLimits(2) yLimits(2)],...
perfColor,'FaceAlpha',0.08,'EdgeColor','none')

yline(PsPermLimit,'--','Color',permColor,'LineWidth',0.8)
yline(PsPerfLimit,'--','Color',perfColor,'LineWidth',0.8)

text(xLimits(1)+0.2,...x
10^( log10(yLimits(1)) + 0.08*(log10(yLimits(2))-log10(yLimits(1))) ),...
'Permeability-limited (E% < 10)','Color',permColor,...
'FontWeight','bold','FontSize',ZONE_FONT)

text(xLimits(1)+0.2,...
10^( log10(yLimits(2)) - 0.08*(log10(yLimits(2))-log10(yLimits(1))) ),...
'Perfusion-limited (E% > 95)','Color',perfColor,...
'FontWeight','bold','FontSize',ZONE_FONT)

text(xLimits(2)-0.2,...
10^( log10(yLimits(2)) - 0.15*(log10(yLimits(2))-log10(yLimits(1))) ),...
'fu_{brp} \times fn_{brp} \times PS_{BBB} = -Q_{cbf} \times ln(1 - E)',...
'HorizontalAlignment','right',...
'FontSize',FORMULA_FONT,'FontWeight','bold','Interpreter','tex')

% -------------------------------------------------------------------------
% Export PDF
% -------------------------------------------------------------------------

set(gcf,'PaperUnits','centimeters')
set(gcf,'PaperPosition',[0 0 FIG_WIDTH FIG_HEIGHT])
set(gcf,'PaperSize',[FIG_WIDTH FIG_HEIGHT])

pdfFile = fullfile(outputDir, 'Figure_3.pdf');
pngFile = fullfile(outputDir, 'Figure_3.png');
print(gcf, pdfFile, '-dpdf', '-vector')
exportgraphics(gcf, pngFile, 'Resolution', 600)

fprintf('Saved %s\n', pdfFile);
fprintf('Saved %s\n', pngFile);

% -------------------------------------------------------------------------
% Local functions
% -------------------------------------------------------------------------

function [c, filled] = mapSourceToStyle(src)

switch char(lower(strtrim(src)))

case 'in situ perfusion'
    c = [0.9255 0.0039 0.0039]; filled = true;

case {'in vitro mdck','ex vivo uc','in vitro caco-2'}
    c = [0 0.4235 0.9451]; filled = true;

case 'in vitro pampa'
    c = [0 0.4745 0]; filled = true;

case 'in silico perfusion'
    c = [0.9255 0.0039 0.0039]; filled = false;

case 'in silico pampa'
    c = [0 0.4745 0]; filled = false;

otherwise
    c = [0 0 0]; filled = true;

end

end
