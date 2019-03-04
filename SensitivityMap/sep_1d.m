function [ind_big,Im]=sep_1d(V,Max_loop)
%This function is used to distinguish two level signal in a vector

%INPUT:
%     V: The image, which is an n X 1 energy vector
%     Max_loop: number of iterations
%OUTPUT
%    BW: the result level set, the value in it is either 1 or -1
%    Im: the mean of image intensity at signal region
%tic

%----------------------------------------------
%    Initialization
%----------------------------------------------
if size(V,1)<size(V,2)
    V = V.';
end
npts = size(V,1);
%Initialize level set
mean_E = mean(V);
phi = ones(npts,1);
phi(find(V<mean_E))=-1;

phi_p1 = find(phi==1);
phi_n1 = find(phi==-1);
m = max(size(phi_p1));
n = npts - m;
H = zeros(npts,1);
H(phi_p1)=1;

Ia = sum(V(phi_p1))/m;
Ib = sum(V(phi_n1))/n;

loop = 0;
%     disp('pause')
%     pause

%----------------------------------------------
%      loop of sweeps
%----------------------------------------------
while loop < Max_loop
    loop = loop +1;
    num_chg = 0;
    %----------------------------------------------
    %    Sweep ROI
    %----------------------------------------------
    for row =1:npts


               
                % -----calculate the differnce of energy-----
                if phi(row)==1
                    
                    dE3 = (V(row) - Ib)^2*n/(n+1) - (V(row) - Ia)^2*m/(m-1);
                elseif phi(row)==-1
                   
                    dE3 = (V(row) - Ia)^2*m/(m+1) - (V(row) - Ib)^2*n/(n-1);
                end  %if phi(row,col)=1
                dE =dE3;
                % -----update according to dE-----
                
                
                if dE<0
                    num_chg = num_chg + 1;
                    if phi(row)==1
                        phi(row) = -1;
                        H(row) = 0;
                        Ia = Ia + (Ia - V(row))/(m-1);
                        Ib = Ib - (Ib - V(row))/(n+1);
                        m = m-1;
                        n = n+1;
                    elseif phi(row)==-1
                        phi(row) = 1;                      
                        H(row) = 1;
                        Ia = Ia - (Ia - V(row))/(m+1);
                        Ib = Ib + (Ib - V(row))/(n-1);
                        m = m+1;
                        n = n-1;
                    end  %if phi(row,col)=1 

                    
                end %if dE<0
                

    end %for row =1:nrow
end%end of while loop < Max_loop
%----------------------------------------------
%      Result and output
%----------------------------------------------
%toc

if Ia<Ib
    ind_big{1} = find(phi==-1);
    Im(1) = Ib;
    ind_big{2} = find(phi==1);
    Im(1) = Ia;
else
    ind_big{1} = find(phi==1);
    Im(1) = Ia;
    ind_big{2} = find(phi==-1);
    Im(2) = Ib;
end
