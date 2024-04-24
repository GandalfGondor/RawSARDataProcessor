function rangeFocusProcessing(imagedata, expParam)

% save image
figure;
imagesc(abs(imagedata(expParam.selectArea(1):expParam.selectArea(2), expParam.selectArea(3):expParam.selectArea(4))));
axis image;colormap('gray');colorbar;
title('Range Compress Image');
ylabel('Range bin');
xlabel('Azimuth bin');
savefig(sprintf('./%s/result/rangeFocus.fig', expParam.fileName));

testData = abs(imagedata(100:1000,2:3));
figure;
imagesc(testData)
axis image
colormap('gray');colorbar;
title('Pulse compression Image')
ylabel('Range bin')
xlabel('Azimuth bin')

% 数据分布
figure;
title('Range Compress Data 幅度分布');
histogram(testData, 100);


% origin
% pdfData = testData(1:200, 1:2000);
pdfData = testData;
pdfData = abs(pdfData(:));
% pdfData = pdfData(1:10000);
% Step 2: 使用 ksdensity 进行核密度估计
[f, x] = ksdensity(pdfData); % f 为估计的概率密度函数值，x 为对应的数据点

% Gamma
gammaPDF = fitdist(pdfData, 'Gamma');
gammaPDFY = pdf(gammaPDF, x);

% Lognormal
lognormalPDF = fitdist(pdfData, 'Lognormal');
lognormalPDFY = pdf(lognormalPDF, x);

% Nakagami
nakagamiPDF = fitdist(pdfData, 'Nakagami');
nakagamiPDFY = pdf(nakagamiPDF, x);

% Rayleigh
rayleighPDF = fitdist(pdfData, 'Rayleigh');
rayleighPDFY = pdf(rayleighPDF, x);

% Weibull
weibullPDF = fitdist(pdfData, 'Weibull');
weibullPDFY = pdf(weibullPDF, x);

% RF-GGD拟合
[theta0, v0, k0] = paramEstimateGGD(pdfData);
lambda = 0.91;
RFGGDPDFY = abs(v0) * k0^k0 / (theta0 * lambda * gamma(k0)) * real((x / theta0) .^ (k0*v0-1)) .* exp(real(-k0 * (x / theta0) .^ v0));

% calculate KL
[gammaNKL,     gammaTKL] = calcKL(f, gammaPDFY);
[lognormalNKL, lognormalTKL] = calcKL(f, lognormalPDFY);
[nakagamiNKL,  nakagamiTKL] = calcKL(f, nakagamiPDFY);
[rayleighNKL,  rayleighTKL] = calcKL(f, rayleighPDFY);
[weibullNKL,   weibullTKL] = calcKL(f, weibullPDFY);
[RFGGDNKL,     RFGGDTKL] = calcKL(f, RFGGDPDFY);

% calculate Kai2
gammaKai     = calcKai(f, gammaPDFY);
lognormalKai = calcKai(f, lognormalPDFY);
nakagamiKai  = calcKai(f, nakagamiPDFY);
rayleighKai  = calcKai(f, rayleighPDFY);
weibullKai   = calcKai(f, weibullPDFY);
RFGGDKai     = calcKai(f, RFGGDPDFY);


% Step 3: 可视化估计的概率密度函数
% 实际数据pdf
figure;
plot(x, f             ,  'ok', 'LineWidth', 1.5, 'Markerfacecolor','k'); hold on;
plot(x, gammaPDFY     , '-*', 'LineWidth', 1.5); hold on;
plot(x, lognormalPDFY ,  '-+', 'LineWidth', 1.5); hold on;
plot(x, nakagamiPDFY  , '-.c', 'LineWidth', 1.5); hold on;
plot(x, rayleighPDFY  , '-s', 'LineWidth', 1.5); hold on;
plot(x, weibullPDFY   , '-v', 'LineWidth', 1.5); hold on;
plot(x, RFGGDPDFY     ,  '-r', 'LineWidth', 1.5);

% set legend
grid on;
xlabel('Normalized Amplitude','FontName','Times NewRoman','FontSize',12);
ylabel('PDF','FontName','Times NewRoman','FontSize',12,'Rotation',0);
lengend1=legend('Histogram', 'Gamma', 'Lognormal', 'Nakagami', 'Rayleigh', 'Weibull', 'RF-GGD', ...
    'Location','Northeast');
set(lengend1,'FontName', 'Times NewRoman', 'FontSize', 12);
axis([0, 800, 0, 5e-3 ]);
savefig(sprintf('./%s/result/linearPdf.fig', expParam.fileName));

% Step 3: 可视化估计的概率密度函数
% 实际数据pdf
figure;
semilogy(x, f             ,  'ok', 'LineWidth', 1.5, 'Markerfacecolor','k'); hold on;
semilogy(x, gammaPDFY     , '-*', 'LineWidth', 1.5); hold on;
semilogy(x, lognormalPDFY ,  '-+', 'LineWidth', 1.5); hold on;
semilogy(x, nakagamiPDFY  , '-.c', 'LineWidth', 1.5); hold on;
semilogy(x, rayleighPDFY  , '-s', 'LineWidth', 1.5); hold on;
semilogy(x, weibullPDFY   , '-v', 'LineWidth', 1.5); hold on;
semilogy(x, RFGGDPDFY     ,  '-r', 'LineWidth', 1.5);

% set legend
grid on;
xlabel('Normalized Amplitude','FontName','Times NewRoman','FontSize',12);
ylabel('PDF','FontName','Times NewRoman','FontSize',12,'Rotation',0);
lengend1=legend('Histogram', 'Gamma', 'Lognormal', 'Nakagami', 'Rayleigh', 'Weibull', 'RF-GGD', ...
    'Location','Southwest');
set(lengend1,'FontName', 'Times NewRoman', 'FontSize', 12);
axis([0, 800, 0, 5e-3 ]);
savefig(sprintf('./%s/result/semilogyPdf.fig', expParam.fileName));

% CFAR
threshold0 = thresholdGGD(lambda, theta0, v0, k0);

CFARGGD(imagedata, threshold0, expParam);


end

