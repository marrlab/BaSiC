lena = im2double(imread('lena.tif'));
figure; imagesc(lena);colorbar; title('original'); % original image
figure; imagesc(mirt_idct2(mirt_dct2(lena))); colorbar; title('reconstruction');% reconstruction from idct2 and dct2


fid = fopen('lena.txt','wt');
for ii = 1:size(lena,1)
    fprintf(fid,'%.8g\t',lena(ii,:));
    fprintf(fid,'\n');
end
fclose(fid)