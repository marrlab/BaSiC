function [A1_hat, E1_hat, A_offset] = inexact_alm_rspca_l1(D, lambda, lambda_darkfield,tol, maxIter, varargin)

% l1 minimizatioin, background has a ratio, rank 1 but not the
% same
% This matlab code implements the inexact augmented Lagrange multiplier 
% method for Sparse low rank matrix recovery

% modified from Robust PCA
% reference:
% Peng et al. "A BaSiC tool for background and shading correction 
% of optical microscopy images" Nature Communications, 14836(2017)
% Candès, E., Li, X., Ma, Y. & Wright, J. "Robust Principal Component
% Analysis?" J. ACM (58) 2011

% D - m x m x n matrix of observations/data (required input)
%

% while ~converged 
%   minimize (inexactly, update A and E only once)
%   L(W, E,Y,u) = |E|_1+lambda * |W|_1 + <Y2,D-repmat(QWQ^T)-E> + +mu/2 * |D-repmat(QWQ^T)-E|_F^2;
%   Y1 = Y1 + \mu * (D - repmat(QWQ^T) - E);
%   \mu = \rho * \mu;
% end
%
% Tingying Peng (tingying.peng@tum.de)

%
% Copyright: CAMP, Technical University of Munich
%       
%% handle varargin
v = varargin;
for i = 1:numel(v)-1
    if ischar(v{i})
            switch lower(v{i})
            case 'weight'
                if numel(v{i+1}) == numel(D)
                    weight = v{i+1};
                else
                    error('weight matrix must be the same size of the image matrix');
                end
            case {'estimatedarkfield'}
                estimatedarkfield = v{i+1}; 
            case{'darkfieldlimit'}
                darkfieldlimit = v{i+1};
            end       
    end
end
  
%%  intialization and given default variables  
p = size(D,1);
q = size(D,2);
m = p*q;
n = size(D,3);
D = reshape(D,m,n);
if exist('weight','var')
    weight = reshape(weight,size(D));
else
    weight = 1;
end
if exist('estimatedarkfield','var')
     if strcmpi(estimatedarkfield,'true')
        darkfield_flag = 1;
    elseif strcmpi(estimatedarkfield,'false')
        darkfield_flag = 0;
    else
        error('"estimatedarkfield" can only be "true" or "false" ');
    end
else
    darkfield_flag = 0;
end 
if exist('darkfieldlimit','var')
    B1_uplimit = darkfieldlimit;
else
    B1_uplimit = 1e7;
end


temp = svd(D, 'econ');
norm_two = temp(1);
clear temp;
Y1 = 0;
Y2 = 0;
ent1 = 1;
ent2 = 10;


A1_hat = zeros(size(D));
E1_hat = zeros(size(D));
W_hat = mirt_dct2(mean(reshape(A1_hat,p,q,n),3));
mu = 12.5/norm_two;% this one can be tuned
mu_bar = mu * 1e7;
rho = 1.5 ;         % this one can be tuned
d_norm = norm(D, 'fro');
A1_coeff = ones(1,n);
A_offset = zeros(m,1);
B1_uplimit = min(D(:));
B1_offset = 0;
A_uplimit = min(D,[],2);
A_inmask = zeros(p,q);
A_inmask(round(p/6):round(p*5/6),round(q/6):round(q*5/6)) = 1;
% main iteration loop starts
iter = 0;
total_svd = 0;
converged = false;
while ~converged       
    iter = iter + 1;
    W_idct_hat = mirt_idct2(W_hat);
    A1_hat = bsxfun(@plus, W_idct_hat(:)*A1_coeff,A_offset);
    temp_W = (D - A1_hat - E1_hat + (1/mu)*Y1 )./ent1;
  
    temp_W = reshape(temp_W,p,q,n);
    temp_W = mean(temp_W,3);
    W_hat = W_hat + mirt_dct2(temp_W);
    W_hat =  max(W_hat - lambda/(ent1*mu),0)+min(W_hat + lambda/(ent1*mu),0);
    W_idct_hat = mirt_idct2(W_hat);
    %W_idct_hat = W_idct_hat;
    A1_hat = bsxfun(@plus, W_idct_hat(:)*A1_coeff,A_offset);
    % update E1 using l0 norm
    E1_hat = E1_hat + (D - A1_hat - E1_hat + (1/mu)*Y1)./ent1;
    E1_hat = max(E1_hat - weight/(ent1*mu), 0)+min(E1_hat + weight/(ent1*mu), 0);
    % update A1_coeff, A2_coeff and A_offset
    %if coeff_flag
        
        %R1 = bsxfun(@minus,D-E1_hat,A_offset); 
        R1 = D-E1_hat;
        A1_coeff = mean(R1)./mean(R1(:));
        %[A1_mincoeff,A1_mincoefftile]= min(A1_coeff);
        %A1_back(A1_coeff<1) = 1;
        A1_coeff(A1_coeff<0) = 0;
        if darkfield_flag == 1
            validA1coeff_idx = find(A1_coeff<1);
            B1_coeff = (mean(R1(W_idct_hat(:)>mean(W_idct_hat(:))-1e-6,validA1coeff_idx))-mean(R1(W_idct_hat(:)<mean(W_idct_hat(:))+1e-6,validA1coeff_idx)))./mean(R1(:));
%             figure(1); plot(A1_coeff);hold on;title('A1_coeff');
%             figure(2); plot(B1_coeff);hold on;title('B1_coeff');
%             temp1 = sum(A1_coeff.^2); temp2 = sum(A1_coeff); 
%             temp3 = sum(B1_coeff); temp4 = sum(A1_coeff.*B1_coeff);
            k = length(validA1coeff_idx);
            temp1 = sum(A1_coeff(validA1coeff_idx).^2); temp2 = sum(A1_coeff(validA1coeff_idx)); 
            temp3 = sum(B1_coeff); temp4 = sum(A1_coeff(validA1coeff_idx).*B1_coeff);
            temp5 = temp2.*temp3-k.*temp4;
            if temp5 == 0
                B1_offset = 0;
            else
                B1_offset = (temp1.*temp3-temp2.*temp4)./temp5;
            end
            % limit B1_offset: 0<B1_offset<B1_uplimit
            B1_offset = max(B1_offset,0);
            B1_offset = min(B1_offset,B1_uplimit./mean(W_idct_hat(:)));
%             %clear temp1 temp2 temp3 temp4;
            B_offset = B1_offset.*mean(W_idct_hat(:))-B1_offset.*W_idct_hat(:);
         %   A_offset = B_offset;
          %   A1_offset = mean(R1,2)-mean(A1_coeff).*W_idct_hat(:);
             A1_offset = mean(R1(:,validA1coeff_idx),2)-mean(A1_coeff(validA1coeff_idx)).*W_idct_hat(:);
             A1_offset = A1_offset-mean(A1_offset(:));
             A_offset = A1_offset-mean(A1_offset(:))-B_offset;
            % smooth A_offset
            W_offset = mirt_dct2(reshape(A_offset,p,q));
            W_offset =  max(W_offset - lambda_darkfield./(ent2*mu),0)+min(W_offset+ lambda_darkfield./(ent2*mu),0);
            A_offset= mirt_idct2(W_offset);
            A_offset = A_offset(:);
            % encourage sparse A_offset
            A_offset = max(A_offset - lambda_darkfield/(ent2*mu), 0)+min(A_offset + lambda_darkfield/(ent2*mu), 0);          
            A_offset = A_offset+B_offset;
        end
        %A_offset= max(min(A_offset,A_uplimit),0);
     Z1 = D - A1_hat - E1_hat;
    %Z2 = D - A2_hat - E2_hat;
    
     Y1 = Y1 + mu*Z1;
    %Y2 = Y2 + mu*Z2; 
    
    mu = min(mu*rho, mu_bar);
        
    %% stop Criterion  
    stopCriterion = norm(Z1, 'fro')/ d_norm;
    %stopCriterion = max(norm(Z1, 'fro')/ d_norm, norm(Z2,'fro')/d_norm);
    if stopCriterion < tol
        converged = true;
    end
    
    if mod( total_svd, 10) == 0
        disp(['Iteration ' num2str(iter) ' |W|_0 ' num2str(sum((abs(W_hat(:))>0)))...
            ' |E1|_0 ' num2str(sum((abs(E1_hat(:))>0)))...
            ' stopCriterion ' num2str(stopCriterion)...
            ' B1_offset ' num2str(B1_offset)]);
    end    

    
%     if mod( total_svd, 10) == 0
%         disp([' |W|_0 ' num2str(sum((abs(W_hat(:))>0)))...
%             ' |E1|_0 ' num2str(sum((abs(E1_hat(:))>0)))...
%             %' |E2|_0 ' num2str(sum((abs(E2_hat(:))>0)))...
%             ' stopCriterion ' num2str(stopCriterion)]);
%     end    
    
    if ~converged && iter >= maxIter
        disp('Maximum iterations reached') ;
        converged = 1 ;       
    end
end
A_offset = A_offset+B1_offset*W_idct_hat(:);

%clear W_hat D weight A_coeff temp_W W_idct_hat Y1 Z1;
