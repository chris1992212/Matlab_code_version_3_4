cd utilities;
fprintf('Compiling by calling mex ...');
mex -O Compute_Ux_Uy.c;
mex -O Compute_Wx_Wy.c;
mex -O Compute_rhs_DxtU_DytU.c;
cd ..;