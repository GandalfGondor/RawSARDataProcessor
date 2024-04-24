function CFARGGD(data, threshold0, expParam)
% CFAR

cfarRangeCompressData = abs(data);
cfar0 = cfarRangeCompressData>threshold0;
% cfar0 = cfar0(selectArea(1):selectArea(2), selectArea(3):selectArea(4));

% % 参数设定
rangeWin = 1;
azimuthWin = 1000;

rStep = 1;
aStep = 200;

[rLen, aLen] = size(cfar0);


record = zeros(rLen, aLen);
for r=1 : rStep : rLen-rangeWin-1
    for a=1 : aStep : aLen-azimuthWin-1
        tmp = cfar0(r: r, a:a+azimuthWin);
        tmp = sum(tmp);
        
        if tmp > expParam.cfarAzimuth % 410
            record(r: r+rangeWin, a:a+azimuthWin) = 1;
        end

        disp(['r=', num2str(r), ', a=', num2str(a)]);
    end
end
figure;
imagesc(record);
axis image;
title('Range Focus CFAR Image');
ylabel('Range bin');
xlabel('Azimuth bin');
savefig(sprintf('./%s/result/CFARGGD.fig', expParam.fileName));

% CFAR效果图
figure;
imagesc(cfar0);
axis image;colormap('gray');colorbar;
title('Range Compress CFAR Result');
ylabel('Range bin');
xlabel('Azimuth bin');
savefig('./result/cfar_result.fig');

end