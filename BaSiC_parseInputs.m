function options = BaSiC_parseInputs(v)
% Parses the command line arguments and fills the options structure with
% default and overridden parameter values 
% 
% Usage:  OPTIONS = parseInputs(v)
%
% Input:  V a cell containing arguments from the BaSiC.m (from varargin)
%
% Output: OPTIONS a structure containing various parameter values needed by
%         BaSiC
%% options that may be specified by the user
options.lambda                  = [];   % default value estimated from input images directly
options.estimation_mode         = [];   % default value = 'l0'
options.max_iterations          = [];   % default value = 500
options.optimization_tol        = [];   % default value = 1e-5;
options.darkfield               = [];   % default value = 'false';
options.lambda_darkfield        = [];   % defulat value estimated from input images directly
% options.basefluo              = [];

%% internal options, should not be reset by user without expert knowledge
options.working_size            = 128;
options.max_reweightiterations     = 10;
options.eplson                  = .1; % reweighting parameter
options.varying_coeff           = 'true';
options.reweight_tol            = 1e-3; % reweighting tolerance


% handle the variable input parameters, provided in (string, param) pairs
for i = 1:numel(v)
    if ischar(v{i})
            switch lower(v{i})
            case 'lambda'
                options.lambda = getParam(v,i);            
            case {'estimation_mode', 'est_mode'}
                options.estimation_mode = getParam(v,i);            
            case 'max_iterations'
                options.max_iterations = getParam(v,i);
            case {'optimization_tol' 'optimisation_tol'}
                options.optimization_tol = getParam(v,i);   
            case {'estimatedarkfield','darkfield'}
                    options.darkfield = getParam(v,i);
            case{'lambda_darkfield','lambda_dark'}
                    options.lambda_darkfield = getParam(v,i);
            case{'basefluo','basefluorescence'}
                    options.basefluo = getParam(v,i);
            end       
    end
    
    if isstruct(v{i})
        options.handles = v{i};
    end
end



% sort the options alphabetically so they are easier to read
options = orderfields(options);



function param = getParam(v,i)

param = [];
if i+1 <= numel(v)
    param = v{i+1};
end