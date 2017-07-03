
function a = mirt_idct2(a)
%  MIRT_IDCT2 2-D  inverse discrete cosine transform; 
%  a is a two dimensional matrix
%  Tingying Peng, 02.03.2015
%  modified from mirt_idctn from Andriy Myronenko

persistent siz ww ind;

%%% Precompute weights and indicies
    siz=size(a);
    ndim=ndims(a);
    if (xor(ndim~=2,~isreal(a))), error('Input requires a real two dimensional matrix');end;
    
    for i=1:ndim,
        n=siz(i); clear tmp;
        
        ww{i} = 2*exp((-1i*pi/(2*n))*(0:n-1)')/sqrt(2*n);
        ww{i}(1) = ww{i}(1)/sqrt(2);
        
        tmp(1:2:n)=(1:ceil(n/2)); 
        tmp(2:2:n)=(n:-1:ceil(n/2)+1);
        ind{i}=bsxfun(@plus, tmp', 0:n:n*(prod(siz)/n-1));
        if (siz(i)==1), break; end;
    end
    



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Actual multidimensional DCT. Handle 1D and 2D cases
%%% separately, because .' is much faster than shiftdim.

% check for 1D or 2D cases

a = idct(idct(a,ww{1},ind{1}).',ww{2},ind{2}).'; % 2D case
    


function a = idct(a,ww,ind)
%DCT  Inverse Discrete cosine transform 1D (operates along first dimension)


    a = bsxfun(@times,  ww, a);  % multiply weights
    a = fft(a);                 % fft
    a=real(a(ind));              % reorder using idicies
    



