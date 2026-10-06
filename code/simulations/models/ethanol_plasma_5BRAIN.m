%ethanol_plasma_5BRAIN plasma + 5BRAIN model of ethanol
%    This function specifies an plasma + 5BRAIN model for predicting
%    ethanol brain concentration.
%    
%    PLASMA PART:
%
%    To execute this model, the following parameters must be defined:
%   
%    Vmax       [Mass/Time]    Maximum velocity of hepatic metabolism
%    Km         [Mass/Volume]  MM constant
%    V1         [Volume]       Central volume of distribution
%    Q          [Mass/Volume]  Intercompartmental clearance
%    V2         [Volume]       Peripheral volume of distribution
%    lambda_po  [1/Time]       First-order oral absorption rate constant
%    F          [unitless]     Oral bioavailability
%    BW_Vper    [unitless]     BW effect on Vper
%
%    BRAIN PART:
%
%    The brain model consists of 4 subspaces:
%       - brain blood (brb); 
%       - brain interstitial space(brint);
%       - brain cellular space (brcel);
%       - cranial cerebrospinal fluid (ccsf);
%       - spinal  cerebrospinal fluid (scsf).
%   
%    To execute this model, the following options must be defined:
%    - tissuePartitioning     
%           a function handle (a tissue partitioning prediction method)


% ========================================================================%
% General structure
% ========================================================================%
function model = ethanol_plasma_5BRAIN()

    % an ODE model definition requires three parts: 
    % initialization, ODEs and output. See below for their definition
    
    model = OdeModel;
    
    model.initfun = @initfun;
    model.rhsfun  = @rhsfun;
    model.obsfun  = @obsfun;
    model.name    = mfilename;
    
end

% ========================================================================%
function setup = initfun(phys, drug, par, options)
%INITFUN Initialization of 12 CMT well-stirred model

    % assertions (model validity):     
    assert(isscalar(drug), 'Model not defined for multiple drugs.')
    

    % indexing
    tissues = {'bra'};
    [queryphys, querydrug, I] = loaddatabases(phys, drug, tissues);
    I = rmfield(I, 'bra');
    I = addcmtidx(I,  'cen', 'per', 'brb', 'GItract', ...
                      'IVbag','IVrate', 'Drinkbag', 'Drinkrate', 'metab', ...
                      'brint','brcel','ccsf','scsf');
    
    Ig = struct;  % groupings
    Ig.organs    = [I.cen I.per I.brb];

    % Define species-specific parameters
    
    % Brain subcompartment volumes supplied from Table S1.
    V_bra         = par.V_brt + par.V_brb;
    V_endothelial = V_bra * par.endothelial_fraction;
    V_csf         = V_bra * par.csf_fraction;
    V.ccsf        = V_csf * par.cranial_csf_fraction;
    V.scsf        = V_csf * (1-par.cranial_csf_fraction);
    V.brb         = par.V_brb;
    V_brm         = V_bra - V_endothelial - V.brb - V.ccsf - V.scsf;
    V.brint       = V_brm * par.brain_ecf_fraction;
    V.brcel       = V_brm * (1-par.brain_ecf_fraction);
    
    % blood flows (ensure closed circulatory system!)
    Q.blo = par.Q_cbf;

    % brain subcompartment flows (according to Verscheijden et al 2019)
    % Maybe set up Ketamine_boost at the script level (via parameters) ?
    Q.productionrate = par.Q_production;
    Q.ssink          = par.spinal_sink_fraction*Q.productionrate;
    Q.csink          = Q.productionrate - Q.ssink;
    Q.sout           = par.spinal_outflow_fraction*Q.ssink;
    Q.sin            = Q.ssink + Q.sout;
    Q.bulk           = par.bulk_flow_fraction*Q.productionrate;
    %Q.met            = par.Qmet;

    % hematocrit
    hct = queryphys('hct');
    
    % Definition of dosing target struct
    Id = struct;
    
    Id.Bolus.iv.cmt     = I.cen;
    Id.Bolus.iv.scaling = 1;       % to be removed later
    
    Id.Infusion.iv.bag  = I.IVbag;
    Id.Infusion.iv.rate = I.IVrate;
    
    Id.Infusion.drink.bag  = I.Drinkbag;
    Id.Infusion.drink.rate = I.Drinkrate;

    Id.Oral.cmt         = I.GItract;
        
    % -------------------------------------------------------------------
    % Define drug-specific parameters

    MW = par.MW;

    % Binding and ionization parameters supplied from Table S3.
    fu.pla   = par.fu_plasma;
    fu.int   = 1;
    fu.cel   = 1;
    fu.brb   = par.fu_brain_plasma/par.BP;
    fu.brint = par.fu_brain_ecf;
    fu.brcel = par.fu_brain_icf;
    fu.ccsf  = par.fu_csf;
    fu.scsf  = par.fu_csf;

    fn.pla   = 1;
    fn.int   = 1;
    fn.cel   = 1;
    fn.brb   = par.fn_brain_plasma;
    fn.brint = par.fn_brain_ecf;
    fn.brcel = par.fn_brain_icf;
    fn.ccsf  = par.fn_csf;
    fn.scsf  = par.fn_csf;
    
    % set brain barriers permeability (change values later)
    PS = struct;
    %PS.b = querydrug('PS_BBB');
    PS.b = par.PS_BBB; 
    PS.c = PS.b * par.BCSFB_to_BBB_ratio;
    PS.e = par.PS_e;
    PS.i = PS.b * par.BCM_to_BBB_ratio;

    % total tissue-to-blood partition coefficients
    BP = par.BP;
    fuB = par.fu_plasma/BP;

    % fraction excreted in feces and metabolized in the gut
    E.gut = querydrug('Egut');
    E.feces = querydrug('Efeces');

    % initial condition and units of ODEs

    doseunit = u.g;
    
    X0 = initializeX0(I);
    X0([I.cen I.per I.brb I.GItract I.IVbag I.Drinkbag I.metab ...
        I.brint I.brcel I.ccsf I.scsf])   = 0 * doseunit;
    X0([I.IVrate I.Drinkrate])  = 0 * doseunit / u.min;
            
    % -----------------------------------------------------------------------
    % Assign model parameters 
    setup = struct;
    setup.indexing.I  = I;
    setup.indexing.Ig = Ig;
    setup.indexing.Id = Id;
    setup.par         = par;
    setup.V           = V;
    setup.Q           = Q;
    setup.fu          = fu;
    setup.fn          = fn;
    setup.PS          = PS;
    setup.BP          = BP;
    setup.MW          = MW;
    setup.hct         = hct;
    setup.fuB         = fuB;
    setup.X0          = X0;
    
    % to be picked up by lump_model 
    setup.cmt         = fieldnames(I);
    [~,setup.physIdx] = ismember(tissues, fieldnames(I));

end



% ========================================================================%
% ODE system of the model
% ========================================================================%
function dX = rhsfun(t, X, setup) % t will be used for infusion rate
%RHSFUN ODE system of 12 CMT well-stirred model
    
    % initialize output vector
    X = X(:);
    dX = NaNsizeof(X);
    
    % model and indexing
    I = setup.indexing.I; 
    %Ig = setup.indexing.Ig;

    % variables (always use column vector notation)
    A_cen     = X(I.cen);
    A_per     = X(I.per);
    A_brb     = X(I.brb);
    A_GItract = X(I.GItract);
    A_brint   = X(I.brint);
    A_brcel   = X(I.brcel);
    A_ccsf    = X(I.ccsf);
    A_scsf    = X(I.scsf);
    infusion_rate = X(I.IVrate);
    drinking_rate = X(I.Drinkrate);

    % tissue volumes, blood flows, extraction ratios, clearance etc.
    V     = setup.V;
    Q     = setup.Q; 

    fuB   = setup.fuB;

    fu = setup.fu;
    fn = setup.fn;
    PS = setup.PS;

    par = setup.par;
    V1  = par.V1;
    V2  = par.V2;
    
    % converting amounts to concentrations
    C_cen  = A_cen / V1;
    C_per  = A_per / V2;
    C_brb  = A_brb / V.brb;
    C_brint = A_brint / V.brint;
    C_brcel = A_brcel / V.brcel;
    C_ccsf = A_ccsf / V.ccsf;
    C_scsf = A_scsf / V.scsf;

    % converting total to unbound and neutral concentration
    Cun_brb   = fn.brb   * fu.brb   * C_brb;
    Cun_brint = fn.brint * fu.brint * C_brint;
    Cun_brcel = fn.brcel * fu.brcel * C_brcel;
    Cun_ccsf  = fn.ccsf  * fu.ccsf  * C_ccsf;
    Cun_scsf  = fn.scsf  * fu.scsf  * C_scsf;



    % ---------------------------------------------------------------------
    % START OF ODEs (in mass units)
    
    % GI tract compartment
    dA_GItract = -par.lambda_po*A_GItract ...
                 + drinking_rate;
        
    % central compartment
    dA_cen        =  par.F * par.lambda_po * A_GItract ...
                   - par.Vmax * C_cen / (par.Km + C_cen) ...
                   + par.Q * (C_per - C_cen) ...
                   + infusion_rate;
    
    % peripheral compartment
    dA_per        =  par.Q * (C_cen - C_per);

    % brain subcompartments
    dA_brb         =   Q.blo*(C_cen    - C_brb) ...
                     + PS.b*(Cun_brint - Cun_brb) ...
                     + PS.c*(Cun_ccsf  - Cun_brb) ...
                     + Q.csink*C_ccsf  + Q.ssink*C_scsf;
    
    dA_brint       = - PS.b*(Cun_brint - Cun_brb) ...
                     - PS.i*(Cun_brint - Cun_brcel) ...
                     - PS.e*(Cun_brint - Cun_ccsf) ...
                     - Q.bulk*C_brint;


    dA_brcel       =   PS.i*(Cun_brint - Cun_brcel); ...

    
    dA_ccsf        =   PS.e*(Cun_brint - Cun_ccsf) ...
                     + PS.c*(Cun_brb - Cun_ccsf) ...
                     + Q.bulk*C_brint ...
                     - Q.sin*C_ccsf + Q.sout*C_scsf ...
                     - Q.csink*C_ccsf;
    
    dA_scsf        =   Q.sin*C_ccsf - Q.sout*C_scsf ...
                     - Q.ssink*C_scsf;
     
    % drug amount in IVbag for infusion
    dA_IVbag   = -infusion_rate;

    % drug amount in Drinkbag for drinking
    dA_Drinkbag = -drinking_rate;
     
    % Change of infusion rate
    d_IVrate   = 0 * unitsOf(X(I.IVrate) / t);

    % Change of drinking rate
    d_Drinkrate   = 0 * unitsOf(X(I.Drinkrate) / t);
     
    % drug amount metabolized or excreted
    dA_metab   = par.Vmax * C_cen / (par.Km + C_cen);
 

    % END OF ODEs 
    % -----------------------------------------------------------------------

    % output vector (always in column vector notation)
    dX(I.cen)        = dA_cen;
    dX(I.per)        = dA_per;
    dX(I.brb)        = dA_brb;
    dX(I.GItract)    = dA_GItract;
    dX(I.metab)      = dA_metab;
    dX(I.IVbag)      = dA_IVbag;
    dX(I.IVrate)     = d_IVrate;
    dX(I.Drinkbag)   = dA_Drinkbag;
    dX(I.Drinkrate)  = d_Drinkrate;
    dX(I.brint)      = dA_brint;
    dX(I.brcel)      = dA_brcel;
    dX(I.ccsf)       = dA_ccsf;
    dX(I.scsf)       = dA_scsf;

end

% ========================================================================%
% Definition of observable quantities
% ========================================================================%
function yobs = obsfun(output, setup, obs)
%OBSFUN Observables in 12 CMT well-stirred model
%   The following observables are supported:
%   
%   Type 'Brain':
%       Site:     'tot','brb','brm','ccsf','sscf'
%       Binding:  'total','unbound'
%       UnitType: 'Mass/Volume','Amount/Volume'
%
%   Type 'PBPK':
%       Site:     'adi','bon','gut','hea','kid','liv','lun','mus','ski','spl','art','ven'
%       Subspace: 'cel','ery','exc','int','pla','tis','tot','vas'
%       Binding:  'total','unbound'
%       UnitType: 'Mass','Amount','Mass/Volume','Amount/Volume'
%   
%   Type 'SimplePK':
%       Site:     'pla'
%       Binding:  'total','unbound'
%       UnitType: 'Mass/Volume','Amount/Volume'
%    
%   Type 'MassBalance':
%       UnitType: 'Mass','Amount'    
    
    I = setup.indexing.I;

    %modelsites = {'adi','bon','gut','hea','kid','liv','lun','mus','ski','spl','art','ven'};
    brainsites = {'tot','brb', 'brint', 'brcel', 'brm','ccsf','scsf'};
    
    yobs = []; % if obs matches a case below, yobs will be overwritten
    
    switch obs.type
        case 'Brain'    % Site, Binding, Unit
            site = obs.attr.Site;
            if ismember(site, brainsites)
                % Site concentration (Mass/Volume units)
                if strcmp(site, 'tot')
                    ytmp = (output.X(:,I.brb) + output.X(:,I.brint) + output.X(:, I.brcel)) / (setup.V.brb + setup.V.brint + setup.V.brcel);
                elseif strcmp(site, 'brm')
                    ytmp = (output.X(:,I.brint) + output.X(:,I.brcel)) / (setup.V.brint + setup.V.brcel);
                else
                    ytmp = output.X(:,I.(site)) / setup.V.(site);
                end
            else
                return
            end
            switch obs.attr.Binding
                case 'total'
                    % do nothing, ytmp already has correct binding type
                case 'unbound'
                    ytmp = ytmp * setup.fu.(site);
                otherwise 
                    return
            end
            switch obs.attr.UnitType 
                case 'Mass/Volume'
                    % do nothing, ytmp already has correct unit type
                case 'Amount/Volume'
                    ytmp = ytmp / setup.MW;
                otherwise 
                    return
            end
            yobs = ytmp;
                
        case 'SimplePK'  % Site, Binding, Unit
            site = obs.attr.Site;
            if strcmp(site, 'pla')
                ytmp = output.X(:,I.cen) / (setup.par.V1 * setup.BP); % bound+unbound, Mass/Volume units
            else
                return
            end
            switch obs.attr.Binding
                case 'total'
                    % do nothing, ytmp already has correct binding type
                case 'unbound'
                    ytmp = ytmp * setup.fu.pla;
                otherwise 
                    return
            end
            switch obs.attr.UnitType 
                case 'Mass/Volume'
                    % do nothing, ytmp already has correct unit type
                case 'Amount/Volume'
                    ytmp = ytmp / setup.MW;
                otherwise 
                    return
            end
            yobs = ytmp;
        case 'SimplePKblood'  % Site, Binding, Unit
            site = obs.attr.Site;
            if strcmp(site, 'blo')
                ytmp = output.X(:,I.cen) / setup.par.V1; % bound+unbound, Mass/Volume units
            else
                return
            end
            switch obs.attr.Binding
                case 'total'
                    % do nothing, ytmp already has correct binding type
                case 'unbound'
                    ytmp = ytmp * setup.fu.pla / setup.BP;
                otherwise 
                    return
            end
            switch obs.attr.UnitType 
                case 'Mass/Volume'
                    % do nothing, ytmp already has correct unit type
                case 'Amount/Volume'
                    ytmp = ytmp / setup.MW;
                otherwise 
                    return
            end
            yobs = ytmp;
 
        case 'MassBalance'
            iIV = I.IVrate;
            ytmp = sum(output.X(:,[1:(iIV-1) (iIV+1):end]),2);
            switch obs.attr.UnitType
                case 'Mass'
                    % do nothing, ytmp already has correct unit type
                case 'Amount'
                    ytmp = ytmp / setup.MW; 
                otherwise
                    return
            end
            yobs = ytmp;
    end
        
end
