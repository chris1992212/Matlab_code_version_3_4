function [ qc ] = QuantCriteria( target, ref, mask )
%QUANTCRI Summary of this function goes here
% qc.mse = MSE4d(target, ref, mask);
qc.rmse = RMSE4d(target, ref, mask);
qc.rr = RelativeError4d(target, ref, mask);
qc.ssim = SSIM4d(target, ref, mask);
end

