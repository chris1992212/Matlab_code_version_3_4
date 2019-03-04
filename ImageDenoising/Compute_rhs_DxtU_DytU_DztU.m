function rhs = Compute_rhs_DxtU_DytU_DztU(Wx,Wy,Wz,bx,by,bz,tau)
rhs = tau*(DxtU(Wx-bx)+DytU(Wy-by) + DztU(Wz-bz));
end

% compute D'_x(U)
function [dxtu] = DxtU(U)
dxtu = [U(:,end,:)-U(:, 1,:) U(:,1:end-1,:)-U(:,2:end,:)];
end

% compute D'_y(U)
function [dytu] = DytU(U)
dytu = [U(end,:,:)-U(1, :,:); U(1:end-1,:,:)-U(2:end,:,:)];
end

function [dztu] = DztU(U)
dztu = zeros(size(U));
dztu(:,:,1) = U(:,:,end) - U(:,:,1);
dztu(:,:,2:end) = U(:,:,1:end-1) - U(:,:,2:end);
end