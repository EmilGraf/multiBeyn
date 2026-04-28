function [points, val, err, cand] = ARMA11(y,opts)
%ARMA11 returns stationary points for the ARMA(1,1) model
%
% [points, val, err, cand] = arma11(y,opts) returns real stationary points
% for the ARMA(1,1) model and objective function ||e||_2^2 
% 
% y(k) + alpha(1)*y(k-1) = e(k) + gamma(1)*e(k-1), k = 2,...,n
%
% for a given vector y of length n.
%
% Input:
%   - y : real vector of size n
%   - opts : options 
%
% Options in opts:
%   - rtol, tolerance for real solutions
%   - restol, residual tolerance
%   - epscluster, tolarance for clustering
%   - xbd,ybd, bounds to search in
%   - zoom, size of contours for first step
%   - p1, quad pts for first step
%
% Output:
%   - points : matrix with real critical points [alpha(1) gamma(1)]
%   - val: values of the objective function 
%   - err: minimal singular values used to verify the solution
%   - cand: matrix with all eigenvalues [alpha(1) gamma(1)]
%
% We find critical values as eigenvalues of a rectangular two-parameter
% eigenvalue problem (A + alpha(1)*B + gamma(1)*C + gamma(2)*D)*z = 0,
% where A,B,C,D are matrices of size (3n-1)*(3n-2)
% Based on arma11 from
% MultiParEig toolbox
% B. Plestenjak and A. Muhic, University of Ljubljana
% P. Holoborodko, Advanpix LLC.
% FreeBSD License, see LICENSE.txt
% Bor Plestenjak 2022
% See: M.E.Hochstenbach, T.Kosir, B.Plestenjak: On the solution of 
% rectangular multiparameter eigenvalue problems, arXiv 2212.01867
% 
% Contour solver from
% multiBeyn
% E. Graf and A. Townsend, Cornell University
% 2026
narginchk(1,2);

% Options
if nargin < 2, opts = []; end
if isfield(opts,'fp_type') && is_numeric_type_supported(opts.fp_type)  
    class_t = opts.fp_type;   
else
    class_t = superiorfloat(y);
end
if ~isfield(opts,'rtol'), opts.rtol = 1e-1;  end
if ~isfield(opts,'restol'), opts.restol = 1e-3;  end
if ~isfield(opts,'xbd'), opts.xbd = [0 1];  end
if ~isfield(opts,'ybd'), opts.ybd = [0 1];  end
if ~isfield(opts,'zoom'), opts.zoom = 2;  end
if ~isfield(opts,'epscluster'), opts.epscluster = 1e-1;  end
if ~isfield(opts,'p1'), opts.p1 = 50;  end

y = y(:);
ytmp = y;
rtol = opts.rtol;
restol = opts.restol;

% we build (3n-1)x(3n-2) matrices such that stationary points are eigenvalues 
% of the rectangular MEP (A + alfa1*B + gamma1*C + gamma1^2*D) 
[A,B,C,D] = ARMA11_matrices(y,class_t);

% use multiBeyn
% sketch
V = randn(size(A,2),size(A,1));
W = randn(size(A,2),size(A,1));
P1 = zeros(size(A,2),size(A,2),3,3);
P2 = zeros(size(A,2),size(A,2),3,3);

P1(:,:,1,1) = V*A;
P1(:,:,2,1) = V*B;
P1(:,:,1,2) = V*C;
P1(:,:,1,3) = V*D;

P2(:,:,1,1) = W*A;
P2(:,:,2,1) = W*B;
P2(:,:,1,2) = W*C;
P2(:,:,1,3) = W*D;

opts.mom = 30;
opts.block = 30;
opts.p = opts.p1;

P121 = P1(:,:,2,1);  P112 = P1(:,:,1,2);  P111 = P1(:,:,1,1);
P122 = P1(:,:,2,2);  P113 = P1(:,:,1,3);  P131 = P1(:,:,3,1);
P221 = P2(:,:,2,1);  P212 = P2(:,:,1,2);  P211 = P2(:,:,1,1);
P222 = P2(:,:,2,2);  P213 = P2(:,:,1,3);  P231 = P2(:,:,3,1);

p1eval = @(z1,z2) P131.*z1.^2 + P121.*z1 + P122.*z1.*z2 + P113.*z2.^2 + P112.*z2 + P111;
p2eval = @(z1,z2) P231.*z1.^2 + P221.*z1 + P222.*z1.*z2 + P213.*z2.^2 + P212.*z2 + P211;

z = opts.zoom;
eigs = [];
for i = opts.xbd(1):2/z:opts.xbd(2)
    for j = opts.ybd(1):2/z:opts.ybd(2)
        x = i; y = j;
        P1(:,:,1,1) = p1eval(x,y);
        P1(:,:,1,2) = (2*P113*y+P122*x+P112)/z;
        P1(:,:,2,1) = (2*P131*x+P122*y+P121)/z;
        P1(:,:,3,1) = P131/z^2;
        P1(:,:,2,2) = P122/z^2;
        P1(:,:,1,3) = P113/z^2;

        P2(:,:,1,1) = p2eval(x,y);
        P2(:,:,1,2) = (2*P213*y+P222*x+P212)/z;
        P2(:,:,2,1) = (2*P231*x+P222*y+P221)/z;
        P2(:,:,3,1) = P231/z^2;
        P2(:,:,2,2) = P222/z^2;
        P2(:,:,1,3) = P213/z^2;
        
        etmp = biBeynQuad(P1,P2,opts);
        if numel(etmp) > 0
            eigs = [eigs; etmp(:,1)/z+x etmp(:,2)/z+y];
        end
    end
end

% filter for real
ind = abs(imag(eigs(:,1))) < rtol & abs(imag(eigs(:,2))) < rtol;
eigs = real(eigs(ind,:));

% --- Clustering ---
einit = cluster_eigs(eigs, opts.epscluster);

% refine
eigs = [];
opts.mom = 20;
opts.block = 20;
opts.p = 50;
z = 40;

for i = 1:size(einit,1)
        x = einit(i,1); y = einit(i,2);
        P1(:,:,1,1) = p1eval(x,y);
        P1(:,:,1,2) = (2*P113*y+P122*x+P112)/z;
        P1(:,:,2,1) = (2*P131*x+P122*y+P121)/z;
        P1(:,:,3,1) = P131/z^2;
        P1(:,:,2,2) = P122/z^2;
        P1(:,:,1,3) = P113/z^2;
        
        P2(:,:,1,1) = p2eval(x,y);
        P2(:,:,1,2) = (2*P213*y+P222*x+P212)/z;
        P2(:,:,2,1) = (2*P231*x+P222*y+P221)/z;
        P2(:,:,3,1) = P231/z^2;
        P2(:,:,2,2) = P222/z^2;
        P2(:,:,1,3) = P213/z^2;

        etmp = biBeynQuad(P1,P2,opts);
        if numel(etmp) > 0
            eigs = [eigs; etmp(:,1)/z+x etmp(:,2)/z+y];
        end
end

% filter for real
ind = abs(imag(eigs(:,1))) < rtol & abs(imag(eigs(:,2))) < rtol;
eigs = real(eigs(ind,:));

% --- Clustering ---
eigs = cluster_eigs(eigs, opts.epscluster);

lambda = eigs;
points = [];
val = [];
err = [];
cand = [];

if numel(lambda) > 0
    alpha = lambda(:,1);
    gamma = lambda(:,2);
    eta = lambda(:,2).^2;
    msvd = zeros(length(alpha),1); 
    for k = 1:length(alpha)
        msvd(k,:) = min(svd(A+alpha(k)*B+gamma(k)*C+eta(k)*D));
    end
    ind = msvd < restol;
    alphaRR = alpha(ind);
    gammaRR = gamma(ind);
    msvdRR = msvd(ind);
    sigmaRR = zeros(length(gammaRR),1);
    y = ytmp(:);
    for k = 1:length(gammaRR)
        sigmaRR(k,1) = arma11_err(y,real(alphaRR(k)),real(gammaRR(k)),class_t);
    end
    points = [alphaRR gammaRR];
    val = sigmaRR;
    cand = [alpha gamma];
    err = msvdRR;
end

end

function [M00,M10,M01,M02] = ARMA11_matrices(y,class_t)
% Returns matrices M00, M10, M01, M11 for the rectangular MEP
%   (M00 + alfa1*M10 + gamma1*M01 + gamma1^2*M02) 
% whose eigenvalues are stationary points of the objective function for the 
% ARMA(1,1) model
 
N = length(y);
R = diag(ones(N-2,1,class_t),1) + diag(ones(N-2,1,class_t),-1);
ZB = zeros(N-1,class_t);
Id = eye(N-1,class_t);
zvec = zeros(N-1,1,class_t);
zrow = zeros(1,N-1,class_t);
y1 = y(1:N-1);
y2 = y(2:N);

M00 = [
    y2    Id    ZB    ZB;
    y1    ZB    Id    ZB;
    zvec  R     ZB    Id; 
    0     y1'   y2'   zrow;
    0     zrow  zrow  y2'
    ];
    
M10 = [
    y1    ZB    ZB    ZB;
    zvec  ZB    ZB    ZB;
    zvec  ZB    ZB    ZB; 
    0     zrow  y1'   zrow;
    0     zrow  zrow  y1'
    ];
    
M01 = [ 
    zvec  R     ZB    ZB;
    zvec  ZB    R     ZB;
    zvec  2*Id  ZB    R;
    0     zrow  zrow  zrow;
    0     zrow  zrow  zrow
    ];
    
M02 = [ 
    zvec  Id    ZB    ZB;
    zvec  ZB    Id    ZB;
    zvec  ZB    ZB    Id;
    0     zrow  zrow  zrow;
    0     zrow  zrow  zrow
    ];
end

function e = arma11_err(y,alpha,gamma,class_t)
% returns value of the objective function for the ARMA(1,1) model
    N = length(y);
    TC = diag(gamma*ones(N-1,1,class_t))+diag(ones(N-2,1,class_t),1); TC(N-1,N) = 1;
    TA = diag(alpha*ones(N-1,1,class_t))+diag(ones(N-2,1,class_t),1); TA(N-1,N) = 1;
    err = pinv(TC)*TA*y;
    e = norm(err)^2;
end

function e_clustered = cluster_eigs(e, tol)
% Helper function to cluster nearby eigenvalues and return their centers
    if isempty(e)
        e_clustered = [];
        return;
    end
    
    n = size(e, 1);
    used = false(n, 1);
    e_clustered = [];
    
    for i = 1:n
        if ~used(i)
            % Compute Euclidean distances for coordinate pairs
            dists = sqrt(sum(abs(e - e(i,:)).^2, 2));
            idx = find(dists < tol);
            
            % Average all points within the threshold tolerance to form cluster center
            e_clustered = [e_clustered; mean(e(idx, :), 1)];
            used(idx) = true;
        end
    end
end