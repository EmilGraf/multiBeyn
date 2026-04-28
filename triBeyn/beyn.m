function [A0,A1] = beyn(M)
% Univariate Beyn
% Returns the Beyn integrals \int f(z) U(P1,P2)
% for f = {1,z2}
% given coefficient tensor M

n = 100; %quadrature nodes;

% Set up quadrature points
quad = 0:2*pi/n:2*pi*(1-1/n);
dquad = exp(2*1i*quad);
quad = exp(1i*quad);

% M is a linear problem
M0 = M{1}; M1 = M{2};
n1 = size(M1,1);

% Calculate matrices for Beyn
A0 = zeros(n1,n1);
A1 = zeros(n1,n1);
for j = 1:n
    A0 = A0 + inv(M0 + M1.*quad(j)).*quad(j);
    A1 = A1 + inv(M0 + M1.*quad(j)).*dquad(j);
end
A0 = A0/n;
A1 = A1/n;

end