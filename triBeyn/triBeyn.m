function [A0,A1] = triBeyn(M)
% Trivariate
% Returns the Beyn integrals \int f(z) U(P1,P2,P3)
% for f = {1,z3}
% given coefficient tensor M

% The contour is the boundary of the unit polydisc:
%   |z1|,|z2|,|z3| ≤ 1.
%
% Each of the three coordinate "faces" (|z_k|=1) is integrated separately.

% Build function handles
P  = cell(1,3);
dP = cell(1,3);

for m = 1:3
    % P_m(z)
    P{m} = @(z) -M{m,1} + M{m,2}*z(1) + M{m,3}*z(2) + M{m,4}*z(3);
    
    % dP_m(z, ell)
    dP{m} = @(z, ell) M{m,ell+1};
end

% Precompute coefficient matrices for Q
Qcoeff = cell(3,4);   % (m, k) where k = 1:4 corresponds to 0, z1, z2, z3

% Get sizes
z0 = zeros(3,1);
P1 = P{1}(z0); P2 = P{2}(z0); P3 = P{3}(z0);

I1 = eye(size(P1,1));
I2 = eye(size(P2,1));
I3 = eye(size(P3,1));

for m = 1:3
    for k = 1:4
        % Extract coefficient matrices from M
        A1 = M{1,k}; A2 = M{2,k}; A3 = M{3,k};
        
        switch m
            case 1
                Qcoeff{m,k} = kron(kron(A1, I2), I3);
            case 2
                Qcoeff{m,k} = kron(kron(I1, A2), I3);
            case 3
                Qcoeff{m,k} = kron(kron(I1, I2), A3);
        end
    end
end

Nr = 20;
Nt = 20;
d = 3;
pref = -2^(d-1)*factorial(d-1)*1i;

% quadrature nodes
[rnodes,w_r] = legpts(Nr,[0,1]);
thetas = (0:Nt-1)*(2*pi/Nt);
w_t = 1/Nt;

% ---- accumulators
A0 = zeros(size(P1,1)*size(P2,1)*size(P3,1));  A1 = A0;

%% ---------- Face k = 1 : |z1| = 1 ----------
for it1 = 1:Nt
    z1 = exp(1i*thetas(it1));
    for it2 = 1:Nt
        th2 = thetas(it2);
        for it3 = 1:Nt
            th3 = thetas(it3);
            for ir2 = 1:Nr
                r2 = rnodes(ir2); w2 = w_r(ir2);
                for ir3 = 1:Nr
                    r3 = rnodes(ir3); w3 = w_r(ir3);
                    z2 = r2*exp(1i*th2);
                    z3 = r3*exp(1i*th3);

                    term = eval_U_tri(1,[z1 z2 z3],Qcoeff);
                    diff_scale = 1i*z1 * r2*r3 * w2*w3 * w_t^3;
                    A0 = A0 + term * diff_scale;
                    A1 = A1 + term * diff_scale * z3;
                end
            end
        end
    end
end

%% ---------- Face k = 2 : |z2| = 1 ----------
for it2 = 1:Nt
    z2 = exp(1i*thetas(it2));
    for it1 = 1:Nt
        th1 = thetas(it1);
        for it3 = 1:Nt
            th3 = thetas(it3);
            for ir1 = 1:Nr
                r1 = rnodes(ir1); w1 = w_r(ir1);
                for ir3 = 1:Nr
                    r3 = rnodes(ir3); w3 = w_r(ir3);
                    z1 = r1*exp(1i*th1);
                    z3 = r3*exp(1i*th3);

                    term = eval_U_tri(2,[z1 z2 z3],Qcoeff);
                    diff_scale = 1i*z2 * r1*r3 * w1*w3 * w_t^3;
                    A0 = A0 - term * diff_scale;
                    A1 = A1 - term * diff_scale * z3;
                end
            end
        end
    end
end

%% ---------- Face k = 3 : |z3| = 1 ----------
for it3 = 1:Nt
    z3 = exp(1i*thetas(it3));
    for it1 = 1:Nt
        th1 = thetas(it1);
        for it2 = 1:Nt
            th2 = thetas(it2);
            for ir1 = 1:Nr
                r1 = rnodes(ir1); w1 = w_r(ir1);
                for ir2 = 1:Nr
                    r2 = rnodes(ir2); w2 = w_r(ir2);
                    z1 = r1*exp(1i*th1);
                    z2 = r2*exp(1i*th2);

                    term = eval_U_tri(3,[z1 z2 z3],Qcoeff);
                    diff_scale = 1i*z3 * r1*r2 * w1*w2 * w_t^3;
                    A0 = A0 + term * diff_scale;
                    A1 = A1 + term * diff_scale * z3;
                end
            end
        end
    end
end

A0 = pref.*A0;
A1 = pref.*A1;
end