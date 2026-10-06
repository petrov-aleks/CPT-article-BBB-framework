%% Generate Figure S13: paracetamol IV versus PO administration

clear; clc;
scriptDir=fileparts(mfilename('fullpath')); repoRoot=fileparts(fileparts(scriptDir));
generatedDir=fullfile(repoRoot,'results','simulations','supplementary_figures','figure_s13');
outputDir=fullfile(repoRoot,'results','figures','supplementary'); if ~exist(outputDir,'dir'), mkdir(outputDir); end
clinical.IV=fullfile(repoRoot,'data','clinical','figure_4','Paracetamol_clinical_data.csv');
clinical.PO=fullfile(repoRoot,'data','clinical','figure_s13','Paracetamol_PO_clinical_data.csv');

FIG_WIDTH=17.8; FIG_HEIGHT=8; AXIS_WIDTH=1.2; SIM_WIDTH=1.5; MARKER_SIZE=5;
DATA_WIDTH=1.2; TICK_FONT=7; LABEL_FONT=7; TITLE_FONT=9; LEGEND_FONT=5.5;
sources=["optimized","perfusion-limited","in situ perfusion","in vitro cell-based", ...
 "in vitro PAMPA","in silico perfusion","in silico PAMPA"];
sites=["plasma","brain ECF","sCSF"]; siteTitles=["Plasma","Brain ECF","CSF"];
colors={[0 0 0],[0 0 0],[.9255 .0039 .0039],[0 .4235 .9451],[0 .4745 0], ...
 [.9255 .0039 .0039],[0 .4745 0]}; styles={'-','--','-','-','-','--','--'};
routes=["IV","PO"];

fig=figure('Units','centimeters','Position',[2 2 FIG_WIDTH FIG_HEIGHT]);
tl=tiledlayout(fig,2,3,'TileSpacing','compact','Padding','compact');
for r=1:2
 route=routes(r); file="Paracetamol_"+route+"_simulations.csv";
 simFile=fullfile(generatedDir,file);
 if ~isfile(simFile)
  error('Missing Figure S13 simulations for %s. Run generate_figure_s13_paracetamol_iv_po_simulations.m first.',route);
 end
 opts=detectImportOptions(simFile); opts=setvartype(opts,{'site','source'},'string'); sim=readtable(simFile,opts);
 sim.source(sim.source=="in vitro MDCK/UC" | sim.source=="in vitro Caco-2")="in vitro cell-based";
 obs=readtable(clinical.(char(route)),'TextType','string');
 for s=1:3
  ax=nexttile(tl); hold(ax,'on'); box(ax,'on'); site=sites(s);
  data=obs(strcmpi(obs.site,site),:); expSources=unique(data.source,'stable');
  markers={'o','^','s'};
  for j=1:min(numel(expSources),3)
   rows=data(data.source==expSources(j),:);
   errorbar(ax,rows.time,rows.conc,rows.sd_down,rows.sd_up,'k','LineStyle','none', ...
    'LineWidth',DATA_WIDTH,'HandleVisibility','off');
   plot(ax,rows.time,rows.conc,markers{j},'MarkerSize',MARKER_SIZE,'MarkerFaceColor','none', ...
    'MarkerEdgeColor','k','LineStyle','none','LineWidth',DATA_WIDTH,'DisplayName',expSources(j));
  end
  rows=sim(strcmpi(sim.site,site),:);
  if site=="plasma"
   plot(ax,rows.time,rows.conc,'k','LineWidth',SIM_WIDTH,'HandleVisibility','off');
   if ~isempty(data), legend(ax,'Location','northeast','FontSize',4,'Box','on'); end
  else
   for j=1:numel(sources)
    sourceRows=rows(rows.source==sources(j),:);
    plot(ax,sourceRows.time,sourceRows.conc,'Color',colors{j},'LineStyle',styles{j}, ...
     'LineWidth',SIM_WIDTH,'HandleVisibility','off');
   end
  end
  set(ax,'LineWidth',AXIS_WIDTH,'FontSize',TICK_FONT); xlim(ax,[0 13]); xticks(ax,0:2:12);
  if route=="IV", plasmaMax=40; brainMax=20; else, plasmaMax=20; brainMax=10; end
  if site=="plasma", ylim(ax,[0 plasmaMax]); else, ylim(ax,[0 brainMax]); end
  xlabel(ax,'Time (h)','FontSize',LABEL_FONT); ylabel(ax,'Concentration (mg/L)','FontSize',LABEL_FONT);
  title(ax,"Paracetamol "+route+" — "+siteTitles(s),'FontSize',TITLE_FONT,'FontWeight','bold');
 end
end

handles=gobjects(numel(sources),1);
for j=1:numel(sources)
 handles(j)=plot(nan,nan,'Color',colors{j},'LineStyle',styles{j},'LineWidth',SIM_WIDTH);
end
lgd=legend(handles,sources,'Orientation','horizontal','FontSize',LEGEND_FONT,'Box','on');
lgd.ItemTokenSize=[20 5]; lgd.Layout.Tile='north';
exportgraphics(fig,fullfile(outputDir,'Figure_S13.pdf'),'ContentType','vector');
exportgraphics(fig,fullfile(outputDir,'Figure_S13.png'),'Resolution',600);
fprintf('Figure S13 saved in %s\n',outputDir);
