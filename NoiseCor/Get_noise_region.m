 function [ind] = Get_noise_region(I,flag) 
 
 %This function is used to define the noise region in an image I
 
 %INPUTs
 %  I: a 2D or 3D matrix
 %  flag: 
 %      1: the inside of the box is noise region
 %      0: the outside of the box is noise region
 %OUTPUT
 %  ind: the index of the noise region
 
 
 if ndims(I)==3
     I = sqrt(sum(I.*conj(I),3));
 end
 
 [nPE nFE] = size(I);
 figure;imagesc(abs(I));axis off;axis equal;colormap(gray)
 
 disp('choose noise box') 
 
 [nxt nyt]=ginput(2);
 
 nx=round(nxt);ny=round(nyt);
 
 Imask = zeros(nPE,nFE);
 Imask(ny(1):ny(2),nx(1):nx(2))=1;
 
 if flag==0
     Imask = 1-Imask;
 end
 
 ind = find(Imask);
 
 close all
 