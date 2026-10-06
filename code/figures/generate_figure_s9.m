%% Generate Figure S9: population concentration-time predictions

clear; clc; close all;

% -------------------------------------------------------------------------
% GLOBAL FIGURE STYLE PARAMETERS
% -------------------------------------------------------------------------

FIG_WIDTH  = 17.8;   % cm
FIG_HEIGHT = 18;     % cm

AXIS_LINE_WIDTH = 1.2;

LINE_WIDTH_MAIN = 1.5;
LINE_WIDTH_DATA = 1.2;

MARKER_SIZE = 5;

X_TICK_FONT = 7;
Y_TICK_FONT = 7;

AXIS_LABEL_FONT = 7;
TITLE_FONT = 9;
PLASMA_LEGEND_FONT = 4;

% -------------------------------------------------------------------------
% Legend location per compound (plasma)
% -------------------------------------------------------------------------

legendLocationMap = containers.Map( ...
    {'Ethanol','Paracetamol','Mannitol','Ibuprofen','Diazepam','Indomethacin'}, ...
    {'northeast','southwest','southeast','southeast','southeast','southeast'} );

% -------------------------------------------------------------------------
% Paths
% -------------------------------------------------------------------------

scriptDir = fileparts(mfilename('fullpath'));
repoRoot = fileparts(fileparts(scriptDir));
clinicalDataDir = fullfile(repoRoot,'data','clinical','figure_4');
generatedSimulationDir = fullfile(repoRoot,'results','simulations', ...
    'supplementary_figures','figure_s9');
outputDir = fullfile(repoRoot,'results','figures','supplementary');
if ~exist(outputDir,'dir'), mkdir(outputDir); end

% -------------------------------------------------------------------------
% Create figure
% -------------------------------------------------------------------------

figure
set(gcf,'Units','centimeters','Position',[2 2 FIG_WIDTH FIG_HEIGHT])

tl = tiledlayout(6,3,'TileSpacing','compact','Padding','compact');

% -------------------------------------------------------------------------
% Compound configuration (same order as main figure)
% -------------------------------------------------------------------------

compounds = {
'Ethanol'      'g/L'     10   1e-2 1e1   3   'brain mass' 0:2:8      [1e-2 1e-1 1e0 1e1]
'Diazepam'     'ng/mL'   25   1e-2 1e3   1   'scsf'       0:4:24     [1e-2 1e0 1e2 1e4 1e6]
'Paracetamol'  'mg/L'    13   1e-2 1e2   3   'scsf'       0:2:12     [1e-2 1e0 1e2]
'Ibuprofen'    'ng/mL'   25   1e-2 1e4   1   'scsf'       0:4:24     [1e-2 1e0 1e2 1e4 1e6]
'Indomethacin' 'ng/mL'   25   1e-2 1e4   1   'scsf'       0:4:24     [1e-2 1e0 1e2 1e4 1e6]
'Mannitol'     'g/L'     7    1e-4 1e1   1   'scsf'       0:2:6      [1e-3 1e-1 1e1]
};

sim_sites = {'plasma','brint','sCSF'};
titles    = {'Plasma','Brain ECF','CSF'};

expMarkers = {'o','^','s'};

tile_id = 1;

% -------------------------------------------------------------------------
% Main loops
% -------------------------------------------------------------------------

for c = 1:size(compounds,1)

    compound      = compounds{c,1};
    unit          = compounds{c,2};
    xmax          = compounds{c,3};
    ymin          = compounds{c,4};
    ymax          = compounds{c,5};
    nExp          = compounds{c,6};
    brainDataSite = compounds{c,7};
    xticks_vals   = compounds{c,8};
    yticks_vals   = compounds{c,9};

    dataFile = fullfile(clinicalDataDir,sprintf('%s_clinical_data.csv',compound));
    simFile = fullfile(generatedSimulationDir,sprintf('%s_population_simulations.csv',compound));
    if ~isfile(dataFile), error('Required clinical data file was not found: %s',dataFile); end
    if ~isfile(simFile)
        error(['Generated population simulation file was not found for %s. ' ...
            'Run generate_figure_s9_population_simulations.m first.'],compound);
    end

    Tdata = readtable(dataFile);
    Tsim  = readtable(simFile);

    for s = 1:3

        ax = nexttile(tile_id);
        ax.Toolbar.Visible = 'off';
        hold(ax,'on')

        sim_name = sim_sites{s};

        % -----------------------------------------------------------------
        % Simulation (no legend)
        % -----------------------------------------------------------------

        Sim = Tsim(strcmpi(Tsim.sim,sim_name),:);

        t   = Sim.time;
        p5  = Sim.p5;
        p50 = Sim.p50;
        p95 = Sim.p95;

        p5(p5<=0)   = 1e-6;
        p50(p50<=0) = 1e-6;
        p95(p95<=0) = 1e-6;

        plot(ax, t, p50,...
            '-k',...
            'LineWidth',LINE_WIDTH_MAIN,...
            'HandleVisibility','off');

        fill(ax,...
            [t; flipud(t)],...
            [p5; flipud(p95)],...
            'k',...
            'FaceAlpha',0.2,...
            'EdgeColor','none',...
            'HandleVisibility','off');

        % -----------------------------------------------------------------
        % Experimental data (goes into legend)
        % -----------------------------------------------------------------

        plotExp = false;

        if strcmp(sim_name,'plasma')
            plotExp  = true;
            dataSite = 'plasma';
        end

        if strcmp(sim_name,'brint') && strcmp(brainDataSite,'brain mass')
            plotExp  = true;
            dataSite = 'brain mass';
        end

        if strcmp(sim_name,'sCSF') && strcmp(brainDataSite,'scsf')
            plotExp  = true;
            dataSite = 'scsf';
        end

        if plotExp

            Tbl        = Tdata(strcmpi(Tdata.site,dataSite),:);
            expSources = unique(Tbl.source);

            for k = 1:min(length(expSources),nExp)

                R = Tbl(strcmp(Tbl.source,expSources{k}),:);

                if ismember('sd_down',Tbl.Properties.VariableNames)

                    errorbar(ax,...
                        R.time,...
                        R.conc,...
                        R.sd_down,...
                        R.sd_up,...
                        'k',...
                        'LineStyle','none',...
                        'LineWidth',LINE_WIDTH_DATA,...
                        'HandleVisibility','off');

                end

                plot(ax,...
                    R.time,...
                    R.conc,...
                    'Marker',expMarkers{k},...
                    'MarkerSize',MARKER_SIZE,...
                    'MarkerFaceColor','none',...
                    'MarkerEdgeColor','k',...
                    'LineStyle','none',...
                    'LineWidth',LINE_WIDTH_DATA,...
                    'DisplayName',expSources{k});

            end

        end

        % -----------------------------------------------------------------
        % Plasma legend ONLY
        % -----------------------------------------------------------------

        if strcmp(sim_name,'plasma')

            lg_loc = legendLocationMap(compound);

            lg = legend(ax,...
                'Location',lg_loc,...
                'FontSize',PLASMA_LEGEND_FONT,...
                'Box','on');

            if strcmp(compound,'Ethanol')
                lg.Position(1) = lg.Position(1)*0.97;
                lg.Position(2) = lg.Position(2)*1.01;
            elseif strcmp(compound,'Paracetamol')
                lg.Position(1) = lg.Position(1)*0.97;
                lg.Position(2) = lg.Position(2)*0.99 ;
            elseif strcmp(compound,'Ibuprofen')
                lg.Position(1) = lg.Position(1)*0.97;
                lg.Position(2) = lg.Position(2)*1.15 ;
            end

        end

        % -----------------------------------------------------------------
        % Axis formatting
        % -----------------------------------------------------------------

        set(ax,'YScale','log','LineWidth',AXIS_LINE_WIDTH)

        ax.XAxis.FontSize = X_TICK_FONT;
        ax.YAxis.FontSize = Y_TICK_FONT;
        ax.YMinorTick = 'off';

        xlim([0 xmax])
        ylim([ymin ymax])

        xticks(ax,xticks_vals)
        yticks(ax,yticks_vals)

        xlabel(ax,'Time (h)','FontSize',AXIS_LABEL_FONT)
        ylabel(ax,['Concentration (' unit ')'],'FontSize',AXIS_LABEL_FONT)

        title(ax,[compound ' — ' titles{s}],...
            'FontSize',TITLE_FONT,'FontWeight','bold')

        box(ax,'on')

        tile_id = tile_id + 1;

    end

end

% -------------------------------------------------------------------------
% Export PDF
% -------------------------------------------------------------------------

set(gcf,'PaperUnits','centimeters')
set(gcf,'PaperPosition',[0 0 FIG_WIDTH FIG_HEIGHT])
set(gcf,'PaperSize',[FIG_WIDTH FIG_HEIGHT])

pdfFile=fullfile(outputDir,'Figure_S9.pdf');
pngFile=fullfile(outputDir,'Figure_S9.png');
print(gcf,pdfFile,'-dpdf','-vector')
exportgraphics(gcf,pngFile,'Resolution',600)
fprintf('Saved %s\n',pdfFile);
fprintf('Saved %s\n',pngFile);
