function [theta0, v0, k0] = paramEstimateGGD(imagedata)

% sampleData = imagedata(1:100,1:100);
sampleData = imagedata;
sampleData = abs(sampleData(:));

sampleDataLen = length(sampleData);

% calc c*
c1 = sum(log(sampleData))/sampleDataLen;

c2 = sum((log(sampleData) - c1).^2)/sampleDataLen;

c3 = sum((log(sampleData) - c1).^3)/sampleDataLen;

% calc a*
a0 = 8*c3^2;

a1 = 4*(3*c3^2-2*c2^3);

a2 = 2*(3*c3^2-8*c2^3);

a3 = c3^2-c2^3;

% q, p
p = (3*a0*a2-a1^2)/(3*a0^2);
q = (2*a1^3 - 9*a0*a1*a2 + 27*a0^2*a3)/(27*a0^3);

% calc K
k0 = (2*c2^3 - 3*c3^2) / (6*c3^2) + ...
    power(-q/2 + sqrt((q/2)^2+(p/3)^3), 1/3) + power(-q/2-sqrt((q/2)^2+(p/3)^3), 1/3);


v0 = sign(-c3)*sqrt((1/c2)*polygamma1(k0));
theta0 = exp(c1- (psi(k0)-log(k0))/v0 );


function res1 = polygamma1(argK)
    res1 = 1 / argK + 1/(2*argK^2);
%     res1 = psi(1, argK);
end

function res2 = polygamma2(argK)
    res2 = -1 / argK^2 - 1/argK^3;
%     res2  = psi(2, argK);
end

end