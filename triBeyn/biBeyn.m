function [A0,A1] = biBeyn(M)
% Bivariate
% Returns the Beyn integrals \int f(z) U(P1,P2)
% for f = {1,z2}
% given coefficient tensor M

% The contour is the boundary of the unit polydisc:
%   |z1|,|z2| ≤ 1.
%
% Each of the 2 coordinate "faces" (|z_k|=1) is integrated separately.

% Build function handles
% Build function handles directly from M
P  = cell(1,3);
dP = cell(1,3);

for m = 1:3
    % P_m(z)
    P{m} = @(z) -M{m,1} + M{m,2}*z(1) + M{m,3}*z(2);
    
    % dP_m(z, ell)
    dP{m} = @(z, ell) M{m,ell+1};
end

% Precompute coefficient matrices for Q
Qcoeff = cell(2,3);   % (m, k) where k = 1:3 corresponds to 0, z1, z2

% Get sizes
z0 = zeros(2,1);
P1 = P{1}(z0); P2 = P{2}(z0);

I1 = eye(size(P1,1));
I2 = eye(size(P2,1));

for m = 1:2
    for k = 1:3
        % Extract coefficient matrices from M
        A1 = M{1,k}; A2 = M{2,k};
        
        switch m
            case 1
                Qcoeff{m,k} = kron(A1, I2);
            case 2
                Qcoeff{m,k} = kron(I1, A2);
        end
    end
end

Nr = 50;
Nt = Nr;
d = 2;
pref = -2^(d-1)*factorial(d-1)*1i;   % (d-1)! / (2π i)^d

% ---- accumulators
A0 = zeros(size(P1,1)*size(P2,1));  A1 = A0;

% ---- quadrature setup
p = Nt;
pts = (1/(2*p)) : (1/p) : (1 - 1/(2*p));
[gpts,gw] = legpts(p,[0,1]);

% precompute exponentials for angles
eAng = exp(2*pi*1i*pts);

% ===================== Integral 1: (z1 on unit circle, z2 radial) =====================
for i = 1:numel(pts)
    z1 = eAng(i);
    for j = 1:numel(pts)
        e2 = eAng(j);
        % (vector over radii)
        for k = 1:numel(gpts)
            r2 = gpts(k);  z2 = r2 * e2;

            term = eval_U_bi(1,[z1 z2],Qcoeff);

            w = gw(k) * ((2/p^3) * z1 * r2);
            A0 = A0 + w * term;                 % f0(z)=1
            A1 = A1 + w * term * z2;
        end
    end
end

% ===================== Integral 2: (z2 on unit circle, z1 radial) =====================
for i = 1:numel(pts)
    e1 = eAng(i);
    for j = 1:numel(pts)
        z2 = eAng(j);
        for k = 1:numel(gpts)
            r1 = gpts(k);  z1 = r1 * e1;

            term = eval_U_bi(2,[z1 z2],Qcoeff);

            w = gw(k) * ((2/p^3) * z2 * r1);
          
            A0 = A0 - w * term;                 % f0(z)=1
            A1 = A1 - w * term * z2;
        end
    end
end

A0 = pref.*A0;
A1 = pref.*A1;
end


