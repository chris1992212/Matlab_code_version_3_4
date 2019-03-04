function [Wx, Wy, Wz] = Compute_Wx_Wy_Wz_PtWs(Ux,Uy,Uz,bx,by,bz,tau)  
%
UUx = Ux + bx; 
UUy = Uy + by;
UUz = Uz + bz;

V = sqrt(UUx.*conj(UUx) + UUy.*conj(UUy) + UUz.*conj(UUz));
V = max(V - tau, 0) ./ max(V,eps);

Wx = V.*UUx; 
Wy = V.*UUy;
Wz = V.*UUz;
