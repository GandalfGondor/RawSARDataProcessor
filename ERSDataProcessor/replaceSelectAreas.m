function tmp = replaceSelectAreas(imagedata, areas)
% 除选定区域均置为zeros

tmp = zeros(size(imagedata));

[x, y] = size(areas);

for i=1:1:x
    tmp(areas(i,1):areas(i,2), areas(i,3):areas(i,4)) = imagedata(areas(i,1):areas(i,2), areas(i,3):areas(i,4));
end


end