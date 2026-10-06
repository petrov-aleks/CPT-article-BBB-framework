%% Generate Figure S1: brain extraction regimes

clear; clc; close all;

% -------------------------------------------------------------------------
% GLOBAL FIGURE STYLE PARAMETERS (CPT:PSP journal compliant)
% -------------------------------------------------------------------------

FIG_WIDTH  = 17.8;    % cm
FIG_HEIGHT = 16;    % cm

AXIS_LINE_WIDTH = 1.2;
CURVE_LINE_WIDTH = 1.2;
CURVE_LINE_WIDTH_INSET = 1;

POINT_LINE_WIDTH = 1.4;

POINT_SIZE = 45;
POINT_SIZE_INSET = 25;

X_TICK_FONT = 10;
Y_TICK_FONT = 10;

AXIS_LABEL_FONT = 10;
TITLE_FONT = 10;
FIG_TITLE_FONT = 12;
LEGEND_FONT = 10;

ZONE_FONT = 7;
INSET_FONT = 6;

% -------------------------------------------------------------------------
% Load data
% -------------------------------------------------------------------------

scriptDir = fileparts(mfilename('fullpath'));
repoRoot = fileparts(fileparts(scriptDir));
tableFile = fullfile(repoRoot, 'results', 'tables', ...
    'figure_3_permeability.csv');
outputDir = fullfile(repoRoot, 'results', 'figures', 'supplementary');

if ~isfile(tableFile)
    error(['Required input was not found: %s\n' ...
        'Run code/tables/generate_table_1.m first.'], tableFile);
end
if ~exist(outputDir, 'dir')
    mkdir(outputDir);
end

T = readtable(tableFile);
T.Compound = string(T.Compound);

% -------------------------------------------------------------------------
% Compounds
% -------------------------------------------------------------------------

compounds = ["Ethanol","Paracetamol","Mannitol", ...
             "Ibuprofen","Diazepam","Indomethacin"];

% -------------------------------------------------------------------------
% fuP and fn parameters
% -------------------------------------------------------------------------

params = table( ...
    compounds', ...
    [1;0.86;1;0.01;0.013;0.01], ...
    [1;1;1;0.0079;1;0.0013], ...
    'VariableNames',{'Compound','fuP','fn'});

% -------------------------------------------------------------------------
% SORT by median(fu × fn × PS) from in situ perfusion
% -------------------------------------------------------------------------

median_eff = zeros(numel(compounds),1);

for i = 1:numel(compounds)

    idx_data   = strcmpi(T.Compound,compounds(i));
    idx_insitu = strcmpi(T.Source,'in situ perfusion');
    idx = idx_data & idx_insitu;

    PS_raw = T.PS_BBB_L_per_h(idx);

    idx_param = strcmpi(params.Compound,compounds(i));
    fu = params.fuP(idx_param);
    fn = params.fn(idx_param);

    PS_eff = fu * fn .* PS_raw;

    median_eff(i) = median(PS_eff,'omitnan');

end

[~,order_idx] = sort(median_eff,'descend');
compounds = compounds(order_idx);

% -------------------------------------------------------------------------
% Brain extraction model
% -------------------------------------------------------------------------

Q = 46.8; % Human cerebral blood flow (L/h)

PS_curve = logspace(-4,6,400);
E_curve  = 100*(1-exp(-PS_curve./Q));

PsPermLimit = 5;
PsPerfLimit = 140;

ExtPermLimit = 100*(1-exp(-PsPermLimit/Q));
ExtPerfLimit = 100*(1-exp(-PsPerfLimit/Q));

permColor = [1 0.5 0];
perfColor = [0.2 0.2 0.2];
alphaZone = 0.08;

% -------------------------------------------------------------------------
% Figure
% -------------------------------------------------------------------------

figure
set(gcf,'Units','centimeters','Position',[2 2 FIG_WIDTH FIG_HEIGHT])

t = tiledlayout(3,2,'Padding','compact','TileSpacing','loose');

title(t,{ ...
    'Brain Extraction (E% = 100 × (1 − e^{-fu_{brp} × fn_{brp} × PS_{BBB}/Q_{CBF}}))'}, ...
    'FontSize',FIG_TITLE_FONT,'FontWeight','bold','Interpreter','tex');

rng(1)

% -------------------------------------------------------------------------
% Loop over compounds
% -------------------------------------------------------------------------

for c = 1:numel(compounds)

    compound = compounds(c);
    ax = nexttile(t);

    % --- Extract data
    idx = strcmpi(T.Compound,compound);
    PS_raw = T.PS_BBB_L_per_h(idx);
    Sources = T.Source(idx);

    idx_param = strcmpi(params.Compound,compound);
    fu = params.fuP(idx_param);
    fn = params.fn(idx_param);

    PS_eff = fu * fn .* PS_raw;
    E_vals = 100*(1-exp(-PS_eff./Q));

    % --- Main Сron–Renkin curve
    semilogx(ax,PS_curve,E_curve,'k-','LineWidth',CURVE_LINE_WIDTH); 
    hold(ax,'on');

    % --- Zones
    xBoxPerm = [1e-6 PsPermLimit PsPermLimit 1e-6];
    yBoxPerm = [0 0 ExtPermLimit ExtPermLimit];

    xBoxPerf = [PsPerfLimit 1e6 1e6 PsPerfLimit];
    yBoxPerf = [ExtPerfLimit ExtPerfLimit 105 105];

    patch(ax,xBoxPerm,yBoxPerm,permColor,...
        'EdgeColor',permColor,'LineStyle',':','FaceAlpha',alphaZone);

    patch(ax,xBoxPerf,yBoxPerf,perfColor,...
        'EdgeColor',perfColor,'LineStyle',':','FaceAlpha',alphaZone);

    % --- Points
    for i = 1:numel(PS_eff)

        [col,filled] = mapSourceToStyle(Sources{i});

        if filled
            scatter(ax,PS_eff(i),E_vals(i),POINT_SIZE,col,'filled');
        else
            scatter(ax,PS_eff(i),E_vals(i),POINT_SIZE,col,...
                'LineWidth',POINT_LINE_WIDTH);
        end

    end

    % --- Axis settings
    set(ax,'XScale','log','LineWidth',AXIS_LINE_WIDTH)

    ax.XAxis.FontSize = X_TICK_FONT;
    ax.YAxis.FontSize = Y_TICK_FONT;
    ax.XMinorTick = 'off';

    xlim(ax,[1e-5 1e6]);
    ylim(ax,[0 105]);
    xticks(ax,10.^(-6:2:8));
    yticks(ax,(0:25:100))
    

    xlabel(ax,'fu_{brp} × fn_{brp} × PS_{BBB} (L/h)','FontSize',AXIS_LABEL_FONT);
    ylabel(ax,'Brain Extraction (%)','FontSize',AXIS_LABEL_FONT);

    title(ax,compound,'FontSize',TITLE_FONT,'FontWeight','bold');

    % --- Zone labels
    text(ax,1e-6*50,ExtPermLimit*1.6,'Permeability-limited',...
        'FontSize',ZONE_FONT,'Color',permColor,'FontWeight','bold');

    text(ax,PsPerfLimit*1.6,ExtPerfLimit*0.93,'Perfusion-limited',...
        'FontSize',ZONE_FONT,'Color',perfColor,'FontWeight','bold');

    % ---------------------------------------------------------------------
    % Inset
    % ---------------------------------------------------------------------
    drawnow

    axPos = ax.Position;

    row = ceil(c/2);

    if row == 1
        y_offset = 0.28;   
    elseif row == 2
        y_offset = 0.32;   
    else
        y_offset = 0.35;   
    end

    inset_pos = [axPos(1) + axPos(3)*0.1, ...
                 axPos(2) + axPos(4)*y_offset, ...
                 axPos(3)*0.35, ...
                 axPos(4)*0.35];



    axInset = axes('Position', inset_pos);

    hold(axInset,'on');
    box(axInset,'on');

    xBoxPerm_in = [1e-6 PsPermLimit PsPermLimit 1e-6];
    yBoxPerm_in = [1e-4 1e-4 ExtPermLimit ExtPermLimit];

    xBoxPerf_in = [PsPerfLimit 1e8 1e8 PsPerfLimit];
    yBoxPerf_in = [ExtPerfLimit ExtPerfLimit 100 100];

    patch(axInset,xBoxPerm_in,yBoxPerm_in,permColor,...
        'EdgeColor',permColor,'LineStyle',':','FaceAlpha',alphaZone);

    patch(axInset,xBoxPerf_in,yBoxPerf_in,perfColor,...
        'EdgeColor',perfColor,'LineStyle',':','FaceAlpha',alphaZone);

    loglog(axInset,PS_curve,E_curve,'k-','LineWidth',CURVE_LINE_WIDTH_INSET);

    for i = 1:numel(PS_eff)

        [col,filled] = mapSourceToStyle(Sources{i});

        if filled
            scatter(axInset,PS_eff(i),E_vals(i),POINT_SIZE_INSET,col,'filled');
        else
            scatter(axInset,PS_eff(i),E_vals(i),POINT_SIZE_INSET,col,...
                'LineWidth',POINT_LINE_WIDTH);
        end

    end

    set(axInset,'XScale','log','YScale','log')
    axInset.XMinorTick = 'off';
    axInset.YMinorTick = 'off';

    xlim(axInset,[1e-4 1e3]);
    ylim(axInset,[1e-4 100]);

    axInset.FontSize = INSET_FONT;
    axInset.XTick = [1e-4 1e-2 1 1e2];
    axInset.YTick = [1e-4 1e-3 1e-2 1e-1 1 10 100];

    title(axInset,'log–log view','FontSize',INSET_FONT);

    hold(axInset,'off');

    set(gcf,'CurrentAxes',ax);

end

% -------------------------------------------------------------------------
% Source styling function
% -------------------------------------------------------------------------

function [c, filled] = mapSourceToStyle(src)

switch lower(strtrim(src))

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

% -------------------------------------------------------------------------
% Global legend
% -------------------------------------------------------------------------

h1 = scatter(NaN,NaN,POINT_SIZE,[0.9255 0.0039 0.0039],'filled');
h2 = scatter(NaN,NaN,POINT_SIZE,[0 0.4235 0.9451],'filled');
h3 = scatter(NaN,NaN,POINT_SIZE,[0 0.4745 0],'filled');
h4 = scatter(NaN,NaN,POINT_SIZE,[0.9255 0.0039 0.0039],'o','LineWidth',POINT_LINE_WIDTH);
h5 = scatter(NaN,NaN,POINT_SIZE,[0 0.4745 0],'o','LineWidth',POINT_LINE_WIDTH);

lgd = legend([h1 h2 h3 h4 h5], ...
{'in situ perfusion','in vitro cell-based','in vitro PAMPA',...
'in silico perfusion','in silico PAMPA'}, ...
'Orientation','horizontal','Box','on','FontSize',LEGEND_FONT);

lgd.Layout.Tile = 'north';

% -------------------------------------------------------------------------
% Export PDF
% -------------------------------------------------------------------------

set(gcf,'PaperUnits','centimeters')
set(gcf,'PaperPosition',[0 0 FIG_WIDTH FIG_HEIGHT])
set(gcf,'PaperSize',[FIG_WIDTH FIG_HEIGHT])

pdfFile = fullfile(outputDir, 'Figure_S1.pdf');
pngFile = fullfile(outputDir, 'Figure_S1.png');
print(gcf, pdfFile, '-dpdf', '-vector')
exportgraphics(gcf, pngFile, 'Resolution', 600)

fprintf('Saved %s\n', pdfFile);
fprintf('Saved %s\n', pngFile);
