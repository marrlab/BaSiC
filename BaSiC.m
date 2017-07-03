function [flatfield, varargout] = BaSiC(IF,varargin)

% Estimation of flatfield for optical microscopy. Apply to a collection of 
% monochromatic images. Multi-channel images should be separated, and each 
% channel corrected separately.
% 
% Usage:  FLATFIELD = BaSiC(IF, ...)
%
% Input:  IF is an NROWSxNCOLZxNIMAG array containing 
%                the images to be corrected.
%
%                BaSiC also supports the following optional arguments,
%                which override default parameters. 
%         
%          'darkfield'(default = 'false') whether you would like to estimate
%          darkfield, keep 'false' if the input images are brightfield
%          images or bright fluoresence images, set 'true' only if only if
%          the input images are fluorescence images are dark and have a strong
%          darkfield contribution. 
%         
%         'basefluo'(default = 'false') set to 'true' if input images has temporal drift (e.g. time lapse movie)      
%
%         'lambda' (default value estimated from input images) high values (eg. 9.5) increase the
%                spatial regularization strength, yielding a more smooth
%                flatfield
%
%          'lambda_dark' (default value estimated from input images) high values (eg. 9.5) increase the
%                spatial regularization strength, yielding a more smooth
%                darkfield
%             
%          'optimization_tol' (default = 1e-5) specifies the tolerance of error in the optimization.
%
%          'max_iterations' (default = 500) specifies the maximum
%                number of iterations allowed in the optimization.
%
% Output: flatfield: estimated flatfield
%         darkfield: estimated darkfield
%         basefluor: estimated per-image background fluoresence signal
% 
% reference: 
% Peng et al. "A BaSiC tool for background and shading correction 
% of optical microscopy images" Nature Communications, 14836(2017) 
%
% March 2016. Tingying Peng: equiry for usage addressed to
% tying84ster@gmail.com or tingying.peng@helmholtz-muenchen.de

addpath('dcttool')

% parse the input arguments, return a structure containing parameters
options = BaSiC_parseInputs(varargin);

% set default values for options that are not specified

if isempty(options.estimation_mode)
    options.estimation_mode = 'l0';
end
if isempty(options.max_iterations)
    options.max_iterations = 500;
end
if isempty(options.optimization_tol)
    options.optimization_tol = 1e-6;
end
if isempty(options.darkfield)
    options.darkfield = 'false';
end
% if isempty(options.basefluo)
%     options.basefluo = 'false';
% end

nrows = options.working_size; ncols = options.working_size; 
D = double(imresize(IF,[nrows ncols],'bilinear'));
%medianD = median(D,3);
meanD = mean(D,3);
meanD = meanD./mean(meanD(:));
W_meanD = mirt_dct2(meanD);

if isempty(options.lambda)
    options.lambda = sum(abs(W_meanD(:)))./(400)*0.5;
end
if isempty(options.lambda_darkfield)
    options.lambda_darkfield = sum(abs(W_meanD(:)))./(400)*0.2;
   % options.lambda_darkfield = max(1e4./sum(abs(minD(:))),0.1);
end


D = sort(D,3);
XAoffset = zeros(nrows,ncols);

weight = ones(size(D));
i = 0;
flag_reweighting = true(1,1);
flatfield_last = ones(nrows,ncols);
darkfield_last = randn(nrows,ncols);
while (flag_reweighting)
       i = i+1;
       disp(['Reweighting Iteration' num2str(i)]); 
       [X_k_A,X_k_E,X_k_Aoffset] = inexact_alm_rspca_l1(D, options.lambda, options.lambda_darkfield,options.optimization_tol,options.max_iterations,'weight',weight,'estimatedarkfield',options.darkfield);
       XA = reshape(X_k_A,nrows,ncols,[]);
       XE = reshape(X_k_E,nrows,ncols,[]);
       XAoffset = reshape(X_k_Aoffset,nrows,ncols);
       XE_norm = XE./(repmat(mean(mean(XA)),nrows,ncols)+1e-6);
       weight = 1./(abs(XE_norm)+options.eplson);
       weight = weight.*numel(weight)./sum(weight(:));
       temp = mean(XA,3)-XAoffset;
       flatfield_current = temp./mean(temp(:));
       darkfield_current = XAoffset;
       mad_flatfield = sum(abs(flatfield_current(:)-flatfield_last(:)))./sum(abs(flatfield_last(:)));
       temp_diff = sum(abs(darkfield_current(:)-darkfield_last(:)));
       if (temp_diff<1e-7)
          mad_darkfield = 0;
       else
          mad_darkfield =temp_diff./(max(sum(abs(darkfield_last(:))),1e-6));
       end
       flatfield_last = flatfield_current;
       darkfield_last = darkfield_current;
       if or(max(mad_flatfield,mad_darkfield)<=options.reweight_tol,i>options.max_reweightiterations)
          flag_reweighting = 0;
       end
end
  

shading =  mean(XA,3)-XAoffset;
%XAoffset = XAoffset+B_offset.*shading;
flatfield = imresize(shading,[size(IF,1) size(IF,2)]);
flatfield = flatfield./mean(flatfield(:));
if strcmpi(options.darkfield,'true')
    darkfield = imresize(reshape(XAoffset,nrows,ncols,[]),[size(IF,1) size(IF,2)]);
else
    darkfield = [];
end
varargout{1} = darkfield;

