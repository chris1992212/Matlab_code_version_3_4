function mask = CreateMask2d(Img)
% 找出头部外延最大区域，生成头部mask
BW = normalize(Img);
% BW = filter2(fspecial('prewitt'),BW);
BW(BW<0.08)=0;
BW(BW>=0.08)=1;

BW = imfill(BW,'holes');

[L,num] = bwlabel(BW);
% 找面积最大的连通区域
if num == 0
    mask = zeros(size(L));
else
    % 找面积最大的连通区域
    TABLE = tabulate(L(:));
    [~,ind] = max(TABLE(2:end,2));
    mask = (L==ind);
end
se = strel('disk',10);
mask = imdilate(mask,se);
mask = imfill(mask,'holes');
mask = imerode(mask,se);

mask = imerode(mask,se);
mask = imdilate(mask,se);
end