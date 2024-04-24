function kai = calcKai(realDist, estiDist)
%CALCKAI 此处显示有关此函数的摘要

 % 剔除非法数据
tmp = find(realDist==0);
tmp = [tmp, find(estiDist==0)];
tmp = unique(tmp,'sorted');
realDist(tmp) = [];
estiDist(tmp) = [];

% 
dataLen = length(realDist);
sumRealDist = sum(realDist);
sumEstiDist = sum(estiDist);

probRealDist = realDist/sumRealDist;
probEstiDist = estiDist/sumEstiDist;

kai = sum((probRealDist - dataLen*probEstiDist).^2 ./ (dataLen*probEstiDist));

end

