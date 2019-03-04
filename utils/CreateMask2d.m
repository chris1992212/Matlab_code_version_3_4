function mask = CreateMask2d(Img)
% �ҳ�ͷ�����������������ͷ��mask
BW = normalize(Img);
% BW = filter2(fspecial('prewitt'),BW);
BW(BW<0.08)=0;
BW(BW>=0.08)=1;

BW = imfill(BW,'holes');

[L,num] = bwlabel(BW);
% �����������ͨ����
if num == 0
    mask = zeros(size(L));
else
    % �����������ͨ����
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