function tableS3 = build_table_s3_indomethacin(baseFile,tableS1File,fitFile,settingsFile)
%BUILD_TABLE_S3_INDOMETHACIN Assemble the Indomethacin column of Table S3.
base=readStrings(baseFile); s1=readtable(tableS1File,'TextType','string');
fit=readtable(fitFile,'TextType','string'); settings=readStrings(settingsFile);
v=@(name) str2double(base.Value(base.Parameter==name)); s1v=@(name) s1.Value(s1.Parameter==name);
sv=@(name) str2double(settings.Value(settings.Parameter==name));
fuCSF=1/(1+v("albumin_csf_to_plasma_ratio")*(1/v("fu_plasma")-1));
fn=@(pH) 1/(1+10^(pH-v("pKa")));
fnBP=fn(v("pH_brain_plasma")); fnECF=fn(v("pH_brain_ecf"));
fnICF=fn(v("pH_brain_icf")); fnCSF=fn(v("pH_csf"));
PSBBB=fit.OptimizedPS_BBB_L_per_h(1); P0=PSBBB/(0.036*s1v("bbb_surface_area"));
PSBCSFB=PSBBB*s1v("bcsfb_to_bbb_surface_area_ratio");
PSBCM=PSBBB*s1v("bcm_to_bbb_surface_area_ratio"); PSE=sv("PS_ependyma");

Category=[repmat("Physicochemical",4,1);repmat("Blood binding",3,1); ...
    repmat("Brain binding",4,1);repmat("Brain ionization",4,1);repmat("Brain permeability",5,1)];
Parameter=["MW";"logP";"Class";"pKa";"BP";"EP";"fu_plasma";"fu_brain_plasma"; ...
    "fu_brain_ECF";"fu_brain_ICF";"fu_CSF";"fn_brain_plasma";"fn_brain_ECF"; ...
    "fn_brain_ICF";"fn_CSF";"P0";"PS_BBB";"PS_BCSFB";"PS_BCM";"PS_E"];
ExactValue=[v("MW");v("logP");NaN;v("pKa");v("BP");v("EP");v("fu_plasma"); ...
    v("fu_brain_plasma");v("fu_brain_ecf");v("fu_brain_icf");fuCSF;fnBP;fnECF;fnICF;fnCSF; ...
    P0;PSBBB;PSBCSFB;PSBCM;PSE];
DisplayValue=[sprintf("%.1f",v("MW"));sprintf("%.2f",v("logP"));"acid";sprintf("%.1f",v("pKa")); ...
    sprintf("%.2f",v("BP"));sprintf("%.3f",v("EP"));sprintf("%.2f",v("fu_plasma")); ...
    sprintf("%.2f",v("fu_brain_plasma"));"1";sprintf("%.3f",v("fu_brain_icf")); ...
    sprintf("%.2f",fuCSF);sprintf("%.4f",fnBP);sprintf("%.4f",fnECF);sprintf("%.4f",fnICF); ...
    sprintf("%.4f",fnCSF);sprintf("%.0f",P0);sprintf("%.0f",PSBBB);sprintf("%.0f",PSBCSFB); ...
    sprintf("%.0f",PSBCM);sprintf("%.0f",PSE)];
Unit=["g/mol";repmat("unitless",14,1);"1e-6 cm/s";repmat("L/h",4,1)];
tableS3=table(Category,Parameter,ExactValue,DisplayValue,Unit);
end

function data=readStrings(file)
opts=detectImportOptions(file); opts=setvartype(opts,opts.VariableNames,'string'); data=readtable(file,opts);
end
