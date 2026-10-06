function obs = FiveBrainobservables(varargin)
%FiveBrainOBSERVABLES Customizable Observable array suitable for brain PBPK modelling.
%   OBS = FiveBrainOBSERVABLES() defines the default brain PBPK observables,
%   a length 9 observable array consisting of plasma concentration and 5(8) 
%   brain sumbcompartments concentrations in Mass/Volume, total concentration.
%
%   OBS = FiveBrainOBSERVABLES(...) allows to customize any of the attributes
%   'Site', 'Binding', 'UnitType' through property-value pairs.
%
%   Examples:
%   
%   obs = FiveBrainobservables();                   % default observables
%   
%   obs = FiveBrainobservables('Binding','Unbound');   % unbound instead of total
%                                              % concentration
%   See also Observable

    % define possible inputs and their defaults
    p = inputParser;

    p.addParameter('Site', {'tot', 'brb', 'brint', 'brcel', 'brm', 'ccsf', 'scsf'})
    p.addParameter('Binding','total');
    p.addParameter('UnitType','Mass/Volume');

    p.parse(varargin{:});

    res = p.Results;

    % create an Observable array from the inputs
    Cpla = Observable('SimplePK','pla',res.Binding,res.UnitType);
    Cbra = Observable('Brain', res.Site, res.Binding, res.UnitType);


    obs = [Cpla; Cbra];

end