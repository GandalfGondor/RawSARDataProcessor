function t0 = thresholdGGD(lambda, theta, v, k)
% 计算 阈值T

t0 = 0;
pfa = 0.08;
if v > 0
    t0 = ((1/k)*gammaincinv(lambda-lambda*pfa, k))^(1/v);
else
    t0 = (1/k)*gammaincinv(lambda*pfa,k)^(1/v);
end

t0 = theta*t0;

end