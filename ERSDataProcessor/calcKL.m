function [normalKL, tailKL] = calcKL(originalRealDist, orginalEstiDist)

% calculate normal KL
normalKL = KLCore(originalRealDist, orginalEstiDist);


% calculate TKL
lenDist = length(originalRealDist);
tailLenDist = round(0.3*lenDist);

tailRealDist = originalRealDist(end-tailLenDist:end);
tailEstiDist = orginalEstiDist(end-tailLenDist:end);
tailKL = KLCore(tailRealDist, tailEstiDist);


% core func
function KL = KLCore(realDist, estiDist)

    % 剔除非法数据
    tmp = find(realDist==0);
    tmp = [tmp, find(estiDist==0)];
    tmp = unique(tmp,'sorted');
    realDist(tmp) = [];
    estiDist(tmp) = [];

    % calculate KL 此处显示有关此函数的摘要
    sumRealDist = sum(realDist);
    sumEstiDist = sum(estiDist);
    
    probRealDist = realDist/sumRealDist;
    probEstiDist = estiDist/sumEstiDist;
    
    factor1 = probRealDist.*log10(probRealDist./probEstiDist);
%     tmp=isinf(factor1) | isnan(factor1);
%     [inf_r inf_c]=find(tmp==1);
%     factor1(inf_c)=[];

    factor2 = probEstiDist.*log10(probEstiDist./probRealDist);
%     tmp=isinf(factor2) | isnan(factor2);
%     [inf_r inf_c]=find(tmp==1);
%     factor2(inf_c)=[];
    
    KL = sum(factor1) + sum(factor2); 
end



end


