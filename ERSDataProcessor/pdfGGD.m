function pdf0 = pdfGGD(lambda, theta0, v0, k0, x)

    lambda = 0.9; 
    x = (0:0.01:1000);
    pdf0 = abs(v0) * k0^k0 / (theta0 * lambda * gamma(k0)) * (x / theta0) .^ (k0*v0-1) .* exp(-k0 * (x / theta0) .^ v0);
    figure;
    plot(x, pdf0);hold on;
    legend;grid on;

end

