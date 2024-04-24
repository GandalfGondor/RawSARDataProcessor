clc;clear all;close all;


%% E1_16687_STD_L0_F93
cfarAzimuth=410;
fileName='E1_16687_STD_L0_F93';
selectArea = [ 873, 895, 11400, 12300;
    110, 130, 9200, 11000;
    430, 440, 8400, 10400;   
    555, 560, 13200, 14200;
];

%%
expParam = ExpParam(selectArea, cfarAzimuth, fileName);
c = physconst('LightSpeed');                % Speed of light

% Extract ERS system parameters
[fs,fc,prf,tau,bw,v,ro,fd] = ERSParameterExtractor(expParam.ldrFileName);

% Extract raw data 
rawData = ERSDataExtractor(expParam.rawFileName,fs,fc,prf,tau,v,ro,fd).';

% Create LFM waveform
waveform = phased.LinearFMWaveform('SampleRate',fs,'PRF',prf,'PulseWidth',tau,'SweepBandwidth',bw,'SweepInterval','Symmetric');
sqang = asind((c*fd)/(2*fc*v));        % Squint angle

% Range migration algorithm
slcimg = myRangeMigrationLFM(rawData,waveform,fc,v,ro,'SquintAngle',sqang, 'ExpParam', expParam);
save('slc.mat', 'slcimg');
Display image
figure;
imagesc(log(abs(slcimg)))
axis image
colormap('gray');colorbar;
title('SLC Image');
ylabel('Range bin');
xlabel('Azimuth bin');


% 数据分布
figure;
title('Range Compress Data 幅度分布');
histogram(abs(slcimg), 100);
% Step 2: 使用 ksdensity 进行核密度估计
[f, xi] = ksdensity(abs(slcimg(:))); % f 为估计的概率密度函数值，xi 为对应的数据点



% save image
slcData = abs(slcimg(selectArea(1):selectArea(2), selectArea(3):selectArea(4)));
imwrite(slcData, sprintf('./result/slc.png'));
figure;
imagesc(abs(slcimg));
axis image
colormap('gray');colorbar;
title('Single Look Complex Image');
ylabel('Range bin');
xlabel('Azimuth bin');
savefig(sprintf('./%s/result/slcTargetFocus.fig', expParam.fileName));
save('slcTargetFocus.mat', 'slcimg');

tmp = zeros(size(slcimg));
tmp(selectArea(1):selectArea(2), selectArea(3):selectArea(4)) = abs(slcimg(selectArea(1):selectArea(2), selectArea(3):selectArea(4)));
slcData = tmp;
figure;
imagesc(slcData);
axis image
colormap('gray');colorbar;
title('Single Look Complex Image');
ylabel('Range bin');
xlabel('Azimuth bin');
savefig('./%s/result/slc.fig');



% Step 3: 可视化估计的概率密度函数
mlimg = multilookProcessing(abs(slcimg),4,20);

% Display Image
figure;
imagesc(mlimg)
axis image
colormap('gray');colorbar;
title('Multi-look Image');
ylabel('Range bin');
xlabel('Azimuth bin');
save('multilook.mat', 'mlimg');