function fi_base = BaSiC_basefluor(IF,flatfield,varargin)

% Estimation of background fluoresence signal for time-lapse movie. Used in
% conjunction with BaSiC.

% Usage:  fi_base = BaSiC(IF, flatfield, ...)
%
% Input:  IF is an NROWSxNCOLZxNIMAG array containing 
%                the images to be corrected.
%
%                BaSiC_basefluor also supports the following optional arguments,
%         
%          'darkfield' to supply your darkfield, note it should be the same size as your fieldfield 

nrows = 128; ncols = 128;
D = double(imresize(IF,[nrows ncols],'bilinear'));
D = reshape(D,nrows*ncols,[]);
flatfield= double(imresize(flatfield,[nrows ncols],'bilinear'));
nVarargs = length(varargin);
if nVarargs ==0
    darkfield = zeros(size(flatfield));
elseif nVarargs == 1
    darkfield = varargin{1};
    darkfield = double(imresize(darkfield,[nrows ncols],'bilinear'));
else
    printf('BaSiC_basefluor accept only one optional input, i.e. darkfield');
end


weight = ones(size(D));
eplson = 0.1;
tol = 1e-6;
for reweighting_iter = 1:5
    W_idct_hat = flatfield(:);
    A_offset = darkfield(:);
    A1_coeff = mean(D);
    % main iteration loop starts
    temp = svd(D, 'econ');
    norm_two = temp(1);
    mu = 12.5/norm_two;% this one can be tuned
    mu_bar = mu * 1e7;
    rho = 1.5 ;         % this one can be tuned
    d_norm = norm(D, 'fro');
    ent1 = 1;
    iter = 0;
    total_svd = 0;
    converged = false;
    A1_hat = zeros(size(D));
    E1_hat = zeros(size(D));
    Y1 = 0;
    while ~converged       
        iter = iter + 1;
        A1_hat = bsxfun(@plus, W_idct_hat*A1_coeff,A_offset);
        % update E1 using l0 norm
        E1_hat = E1_hat + (D - A1_hat - E1_hat + (1/mu)*Y1)./ent1;
        E1_hat = max(E1_hat - weight/(ent1*mu), 0)+min(E1_hat + weight/(ent1*mu), 0);
        % update A1_coeff, A2_coeff and A_offset
        %if coeff_flag

            %R1 = bsxfun(@minus,D-E1_hat,A_offset); 
            R1 = D-E1_hat;
            A1_coeff = mean(R1)-mean(A_offset);
            %[A1_mincoeff,A1_mincoefftile]= min(A1_coeff);
            %A1_back(A1_coeff<1) = 1;
            A1_coeff(A1_coeff<0) = 0;

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
            disp(['Iteration ' num2str(iter) ...
                ' |E1|_0 ' num2str(sum((abs(E1_hat(:))>0)))...
                ' stopCriterion ' num2str(stopCriterion)]);
        end  
    end
    % update weight
    XE_norm = bsxfun(@ldivide, E1_hat, mean(A1_hat));
    weight = 1./(abs(XE_norm)+eplson);
%     if isempty(segmentation)
%         weight(segmentation~=0)=0;
%     end
    weight = weight.*numel(weight)./sum(weight(:));    
end
fi_base = A1_coeff;
