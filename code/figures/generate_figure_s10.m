%% Generate Figure S10: influence of CNS surface-area assumptions

clear; clc;
scriptDir=fileparts(mfilename('fullpath'));
repoRoot=fileparts(fileparts(scriptDir));
generatedDir=fullfile(repoRoot,'results','simulations','supplementary_figures','figure_s10');
scenarioFile=fullfile(repoRoot,'data','model_parameters','figure_s10_surface_area_scenarios.csv');
outputDir=fullfile(repoRoot,'results','figures','supplementary');
if ~exist(outputDir,'dir'), mkdir(outputDir); end

FIG_WIDTH=17.8; FIG_HEIGHT=18; AXIS_LINE_WIDTH=1.2; LINE_WIDTH_SIM=1.5;
X_TICK_FONT=7; Y_TICK_FONT=7; AXIS_LABEL_FONT=7; TITLE_FONT=9; LEGEND_FONT=5.5;
config={
 'Ethanol','g/L',9,1e-2,1e1,0:2:8,[1e-2 1e-1 1e0 1e1]
 'Diazepam','ng/mL',25,1e-2,1e3,0:4:24,[1e-2 1e0 1e2 1e4 1e6]
 'Paracetamol','mg/L',13,1e-2,1e2,0:2:12,[1e-2 1e0 1e2]
 'Ibuprofen','ng/mL',25,1e-2,1e4,0:4:24,[1e-2 1e0 1e2 1e4 1e6]
 'Indomethacin','ng/mL',25,1e-2,1e4,0:4:24,[1e-2 1e0 1e2 1e4 1e6]
 'Mannitol','g/L',7,1e-4,1e1,0:2:6,[1e-3 1e-1 1e1]};
sites=["brain ICF","brain ECF","sCSF"]; titles=["Brain ICF","Brain ECF","CSF"];
scenarios=readtable(scenarioFile,'TextType','string'); colors=lines(height(scenarios));

fig=figure('Units','centimeters','Position',[2 2 FIG_WIDTH FIG_HEIGHT]);
tl=tiledlayout(fig,6,3,'TileSpacing','compact','Padding','compact'); tileId=1;
for c=1:size(config,1)
 compound=string(config{c,1}); file=compound+"_surface_area_simulations.csv";
 inputFile=fullfile(generatedDir,file);
 if ~isfile(inputFile)
  error('Missing Figure S10 simulations for %s. Run generate_figure_s10_surface_area_simulations.m first.',compound);
 end
 opts=detectImportOptions(inputFile); opts=setvartype(opts,{'site','PS'},'string'); sim=readtable(inputFile,opts);
 for s=1:3
  ax=nexttile(tl,tileId); hold(ax,'on'); siteData=sim(strcmpi(sim.site,sites(s)),:);
  for j=1:height(scenarios)
   rows=siteData(siteData.PS==scenarios.Scenario(j),:);
   plot(ax,rows.time,rows.conc,'Color',colors(j,:),'LineWidth',LINE_WIDTH_SIM,'HandleVisibility','off');
  end
  set(ax,'YScale','log','LineWidth',AXIS_LINE_WIDTH); ax.XAxis.FontSize=X_TICK_FONT;
  ax.YAxis.FontSize=Y_TICK_FONT; ax.YMinorTick='off';
  xlim(ax,[0 config{c,3}]); ylim(ax,[config{c,4} config{c,5}]);
  xticks(ax,config{c,6}); yticks(ax,config{c,7});
  xlabel(ax,'Time (h)','FontSize',AXIS_LABEL_FONT);
  ylabel(ax,"Concentration ("+string(config{c,2})+")",'FontSize',AXIS_LABEL_FONT);
  title(ax,compound+" — "+titles(s),'FontSize',TITLE_FONT,'FontWeight','bold'); box(ax,'on'); tileId=tileId+1;
 end
end

handles=gobjects(height(scenarios),1); labels=strings(height(scenarios),1);
for j=1:height(scenarios)
 handles(j)=plot(nan,nan,'Color',colors(j,:),'LineWidth',LINE_WIDTH_SIM);
 labels(j)=sprintf('SA_{BBB} = %g m^2, SA_{BCSFB} = %g m^2, SA_{BCM} = %g m^2', ...
  scenarios.SA_BBB_m2(j),scenarios.SA_BCSFB_m2(j),scenarios.SA_BCM_m2(j));
end
lgd=legend(handles,labels,'Orientation','horizontal','FontSize',LEGEND_FONT,'Box','on');
lgd.Layout.Tile='north'; lgd.NumColumns=2;
exportgraphics(fig,fullfile(outputDir,'Figure_S10.pdf'),'ContentType','vector');
exportgraphics(fig,fullfile(outputDir,'Figure_S10.png'),'Resolution',600);
fprintf('Figure S10 saved in %s\n',outputDir);
