function eigs = biBeynQuad(P1,P2,opts)
% Bivariate Beyn's method based on Grothendieck residue, quadratic version

% Usage: eigs = biBeyn2(P1, P2, opts)
% - P1, P2 are the 3D coefficient tensors for a quadratic 2PEP
% - opts (optional): struct with fields
% .tol (default 1e-4)
% .p (Gauss-Legendre nodes per radial; default 50)
% .block (random test block size; default k at the moment)
% .h (height for ellipse

if nargin < 3, opts = struct; end
if ~isfield(opts,'tol'), opts.tol = 1e-4; end
if ~isfield(opts,'p'), opts.p = 50; end
if ~isfield(opts,'mom'), opts.mom = 5; end
if ~isfield(opts,'h'), opts.h = 1e-6; end

% dimensions
n1 = size(P1,1); n2 = size(P2,1); k = n1*n2;
if ~isfield(opts,'block'), opts.block = ceil(min(k,64)/opts.mom); end

u1 = randn(1,n1,opts.block) + 1i*randn(1,n1,opts.block);
v1 = randn(n2,1,opts.block) + 1i*randn(n2,1,opts.block);
u2 = randn(1,n1,opts.block) + 1i*randn(1,n1,opts.block);
v2 = randn(1,n2,opts.block) + 1i*randn(1,n2,opts.block);

P121 = P1(:,:,2,1);  P112 = P1(:,:,1,2);  P111 = P1(:,:,1,1);
P122 = P1(:,:,2,2);  P113 = P1(:,:,1,3);  P131 = P1(:,:,3,1);
P221 = P2(:,:,2,1);  P212 = P2(:,:,1,2);  P211 = P2(:,:,1,1);
P222 = P2(:,:,2,2);  P213 = P2(:,:,1,3);  P231 = P2(:,:,3,1);

[A0,A1,A2] = beyn(P111,P112,P121,P122,P113,P131,P211,P212,P221,P222,P213,P231,u1,v1,u2,v2,opts);

eigs = matcheigs(A0,A1,A2,P111,P112,P121,P122,P113,P131,P211,P212,P221,P222,P213,P231,opts,1e-1);

end

% =========================================================================
function [A,B,C] = beyn(P111,P112,P121,P122,P113,P131,P211,P212,P221,P222,P213,P231,u1,v1,u2,v2,opts)
% Returns the Beyn integrals \int W f(z) U(P1,P2) V
% for f in {1,z1^ii,z2ii}, ii = 1,...,2*opts.mom-1

% ---- quadrature setup (mapping of poly-ellipse boundary to [0,1]^3)
p = opts.p;
m = 2*opts.mom-1;
b = opts.h;
pts = (1/(2*p)) : (1/p) : (1 - 1/(2*p));
[gpts,gw] = legpts(p,[0,1]);
ew = 1i*b*cos(2*pi*pts) - sin(2*pi*pts);

% precompute trig for angles
Cos = cos(2*pi*pts);
Sin = b*sin(2*pi*pts);

% identities (reused inside kronapp calls)
I1 = eye(size(P111,1));
I2 = eye(size(P211,1));

% ---- accumulators
s = [size(u1,3) size(u1,3) m];
A = zeros(s(1:2));  B = zeros(s);  C = B;

% ---- Jacobian partials
d1z1 = @(z1,z2) 2*P131*z1 + P121 + P122*z2;
d2z1 = @(z1,z2) 2*P231*z1 + P221 + P222*z2;
d1z2 = @(z1,z2) 2*P113*z2 + P112 + P122*z1;
d2z2 = @(z1,z2) 2*P213*z2 + P212 + P222*z1;


d1z1 = @(z1,z2,V) kronapp(d1z1(z1,z2)',eye(size(P211,1)),V);
d1z2 = @(z1,z2,V) kronapp(d1z2(z1,z2)',eye(size(P211,1)),V);
d2z1 = @(z1,z2,V) kronapp(eye(size(P111,1)),d2z1(z1,z2)',V);
d2z2 = @(z1,z2,V) kronapp(eye(size(P111,1)),d2z2(z1,z2)',V);

% Evaluate P1,P2 quickly at (z1,z2)
p1eval = @(z1,z2) P131*z1^2 + P121*z1 + P122*z1*z2 + P113*z2^2 + P112*z2 + P111;
p2eval = @(z1,z2) P231*z1^2 + P221*z1 + P222*z1*z2 + P213*z2^2 + P212*z2 + P211;

% ===================== Integral 1: (z1 on unit circle, z2 radial) =====================
parfor i = 1:numel(pts)
    Bi = zeros(s);
    Ci = zeros(s);
    z1 = Cos(i)+1i*Sin(i);
    for j = 1:numel(pts)
        e2 = Cos(j)+1i*Sin(j);
        % (vector over radii)
        for k = 1:numel(gpts)
            r2 = gpts(k);  z2 = r2 * e2;

            % evaluate once, reuse
            P1z = p1eval(z1,z2);
            P2z = p2eval(z1,z2);

            % core Lyapunov/Sylvester pieces (reuse factors)
            Qa  = P1z*P1z';
            Qd  = P2z*P2z';
            Q1a = P1z'*P1z;
            % Q1d = Qd;                % = P2z*P2z'
            % Q2a = Qa;                % = P1z*P1z'
            Q2d = P2z'*P2z;

            % precalculate Schur factors
            [QaZ,QaT]  = schur(Qa.','complex');
            [QdZ,QdT]  = schur(Qd,'complex');
            [Q1aZ,Q1aT]  = schur(Q1a,'complex');
            [Q1dZ,Q1dT]  = schur(Qd.','complex');
            [Q2aZ,Q2aT]  = schur(Qa,'complex');
            [Q2dZ,Q2dT]  = schur(Q2d.','complex');

            J = syl(QaZ,QaT,QdZ,QdT,u1,v1);

            % J11 = (P1z')* syl(Q2a,Q2d, d2z2(J))  (kronapp with (P1z') ⊗ I2)
            % J11 = kronapp(P1z', I2, syl(Q2aZ,Q2aT,Q2dZ,Q2dT,d2z2(J)));
            [x1,x2] = kronappk(P1z',I2,u2,v2);
            
            J11 = syl(Q2aZ,Q2aT,Q2dZ,Q2dT,x1,pagetranspose(x2)).'*d2z2(z1,z2,J);

            % J12 = (I1 ⊗ P2z') * syl(Q1a,Q1d, d1z2(J))
            % J12 = kronapp(I1, P2z', syl(Q1aZ,Q1aT,Q1dZ,Q1dT,d1z2(J)));
            [x1,x2] = kronappk(I1,P2z',u2,v2);

            J12 = syl(Q1aZ,Q1aT,Q1dZ,Q1dT,x1,pagetranspose(x2)).'*d1z2(z1,z2,J);

            J = J11 - J12;

            w = gw(k) * (ew(i) * r2);
            A = A + w * J;                 % f0(z)=1
            for ii = 1:m
                Bi(:,:,ii) = Bi(:,:,ii) + w * (z1^ii * J);          % f1(ii)(z)=z1^ii
                Ci(:,:,ii) = Ci(:,:,ii) + w * (z2^ii * J);          % f2(ii)(z)=z2^ii
            end
        end
    end
    B = B + Bi;
    C = C + Ci;
end

% ===================== Integral 2: (z2 on unit circle, z1 radial) =====================
parfor i = 1:numel(pts)
    Bi = zeros(s);
    Ci = zeros(s);
    e1 = Cos(i)+1i*Sin(i);
    for j = 1:numel(pts)
        z2 = Cos(j)+1i*Sin(j);
        for k = 1:numel(gpts)
            r1 = gpts(k);  z1 = r1 * e1;

            P1z = p1eval(z1,z2);
            P2z = p2eval(z1,z2);

            Qa  = P1z*P1z';
            Qd  = P2z*P2z';
            Q1a = P1z'*P1z;
            % Q1d = Qd;                % = P2z*P2z'
            % Q2a = Qa;                % = P1z*P1z'
            Q2d = P2z'*P2z;
            
            % precalculate Schur factors
            [QaZ,QaT]  = schur(Qa.','complex');
            [QdZ,QdT]  = schur(Qd,'complex');
            [Q1aZ,Q1aT]  = schur(Q1a,'complex');
            [Q1dZ,Q1dT]  = schur(Qd.','complex');
            [Q2aZ,Q2aT]  = schur(Qa,'complex');
            [Q2dZ,Q2dT]  = schur(Q2d.','complex');

            J = syl(QaZ,QaT,QdZ,QdT,u1,v1);

            % J21 = (I1 ⊗ P2z') * syl(Q1a,Q1d, d1z1(J))
            % J21 = kronapp(I1, P2z', syl(Q1aZ,Q1aT,Q1dZ,Q1dT,d1z1(J)));
            [x1,x2] = kronappk(I1,P2z',u2,v2);
            
            J21 = syl(Q1aZ,Q1aT,Q1dZ,Q1dT,x1,pagetranspose(x2)).'*d1z1(z1,z2,J);

            % J22 = (P1z') ⊗ I2 * syl(Q2a,Q2d, d2z1(J))
            % J22 = kronapp(P1z', I2, syl(Q2aZ,Q2aT,Q2dZ,Q2dT,d2z1(J)));  
            [x1,x2] = kronappk(P1z',I2,u2,v2);

            J22 = syl(Q2aZ,Q2aT,Q2dZ,Q2dT,x1,pagetranspose(x2)).'*d2z1(z1,z2,J);

            J = J21 - J22;

            w = gw(k) * (ew(j) * r1);
            A = A + w * J;
            for ii = 1:m
                Bi(:,:,ii) = Bi(:,:,ii) + w * (z1^ii * J);          % f1(ii)(z)=z1^ii
                Ci(:,:,ii) = Ci(:,:,ii) + w * (z2^ii * J);          % f2(ii)(z)=z2^ii
            end
        end
    end
    B = B + Bi;
    C = C + Ci;
end
end

% =========================================================================
% Fast kron-apply: kron(A,B)*C with many-right-hand-sides in C
function X = kronapp(A,B,C)
m = size(B,1); n = size(A,1);
X = reshape(C,m,n,[]);
X = pagemtimes(pagemtimes(B,X),A.');   % batched B*X*A.'
X = reshape(X,n*m,[]);
end

% kron-apply C*kron(A,B) where rows of C = Cu⊗Cv have Kronecker structure
% keep structure in the result
function [x1,x2] = kronappk(A,B,Cu,Cv)
x1 = pagemtimes(Cu,A);
x2 = pagemtimes(Cv,B);
end

%%% =========================================================================
% Fast Sylvester with precomputed Schur for positive definite A,B
% Schur form is always diagonal
% Assuming Kronecker structured C
function X = syl(ZB,TB,ZA,TA,Cu,Cv)
% Solves (I ⊗ A + B ⊗ I) \ C  with multiple RHS (columns of C)
% C(:,:,i) = Cu(:,:,i) ⊗ Cv(:,:,i)

m = size(ZA,1); n = size(ZB,1); k = size(Cu,3);

TA_diag = diag(TA); TB_diag = diag(TB);

% -- Transform RHS: F = ZA' * C * ZB
F = pagemtimes(pagemtimes(ZA,'ctranspose',Cv,'none'),pagemtimes(Cu,ZB));

% =========================
% FAST PATH: both diagonal
% =========================
% Y_ij,: = -F_ij,: ./ (dA_i + dB_j)
denom = TA_diag + TB_diag.';             % m x n
Y = F ./ denom;                          % implicit expansion to m x n x k


% -- Back transform
X = reshape(pagemtimes(ZA,pagemtimes(Y,'none',ZB,'ctranspose')), m*n, k);
end

%%% =========================================================================
% Hankel Eigensolver and plug in and solve again
function eigs = matcheigs(A0,A1,A2,P111,P112,P121,P122,P113,P131,P211,P212,P221,P222,P213,P231,opts,tol)
    % build Hankel matrices
    m = opts.mom;
    ind = size(A0,1);
    H11 = zeros(m*ind); H10 = H11; H21 = H11; H20 = H11;
    for i = 0:m-1
        for j = 0:m-1
            H11(1+i*ind:(i+1)*ind,1+j*ind:(j+1)*ind) = A1(:,:,i+j+1);
            H21(1+i*ind:(i+1)*ind,1+j*ind:(j+1)*ind) = A2(:,:,i+j+1);
            if i+j > 0
                H10(1+i*ind:(i+1)*ind,1+j*ind:(j+1)*ind) = A1(:,:,i+j);
                H20(1+i*ind:(i+1)*ind,1+j*ind:(j+1)*ind) = A2(:,:,i+j);
            else
                H10(1+i*ind:(i+1)*ind,1+j*ind:(j+1)*ind) = A0;
                H20(1+i*ind:(i+1)*ind,1+j*ind:(j+1)*ind) = A0;
            end
        end
    end
    
    % project
    [U,S,Vsvd] = svd(H10,'econ');
    s = diag(S);
    ind = nnz(s > opts.tol*s(1));
    if ind == 0, ind = min(10,size(S,2)); end % fall back modest rank
    H11 = U(:,1:ind)'*H11*Vsvd(:,1:ind);
    H10 = S(1:ind,1:ind);
    
    % get eigenvalues and eigenvectors
    e1 = eig(H11,H10);
    eigs = [];

    p1eval = @(z1,z2) P131*z1^2 + P121*z1 + P122*z1*z2 + P113*z2^2 + P112*z2 + P111;
    p2eval = @(z1,z2) P231*z1^2 + P221*z1 + P222*z1*z2 + P213*z2^2 + P212*z2 + P211;

    % plug in and solve subproblem
    for i = 1:size(e1)
        e21 = polyeig(P131*e1(i).^2+P121*e1(i)+P111,P122*e1(i)+P112,P113);
        e22 = polyeig(P231*e1(i).^2+P221*e1(i)+P211,P222*e1(i)+P212,P213);
        i1 = abs(e21) < 1;
        i2 = abs(e22) < 1;
        e21 = e21(i1);
        e22 = e22(i2);
        i1 = zeros(size(e21));
        if numel(e21) > 0 && numel(e22) > 0
            for j = 1:size(e21)
                [i1(j),ind] = min(abs(e21(j)-e22),[],"all");
                e21(j) = (e21(j) + e22(ind))/2;
            end
            i1 = (i1 < tol);
            e2 = e21(i1);
            % residual check
            eval = zeros(numel(e2),1);
            for j = 1:numel(e2)
                eval(j) = min(svd(p1eval(e1(i),e2(j))))+min(svd(p2eval(e1(i),e2(j))));
            end
            [val,ind] = min(eval);
            if val < tol
                eigs = [eigs; e1(i) e2(ind)];
            end
        end
    end

    %%% trim outside poly-ellipse
    if size(eigs,1) > 0
        z1 = eigs(:,1);
        z2 = eigs(:,2);

        ind = ( (real(z1).^2 + (imag(z1)).^2) < 1 ) & ...
            ( (real(z2).^2 + (imag(z2)).^2) < 1 );

        eigs = (eigs(ind,:));
    end
end




