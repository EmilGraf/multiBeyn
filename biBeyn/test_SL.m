%%%Test multiparameter Sturm-Liouville with asymptotics
% Asymptotics
asy = [
    13498.430479   4642.853415;
    18028.248411   13606.167147;
    22577.169734   26842.168609;
    53993.721914   18571.413659;
    72112.993645   54424.668586
];
eigenpairs = [];

% find eigenpairs using biBeyn
for i = 1:5
    rng(1);
    eigenpairs = [eigenpairs; testmpsl(asy(i,1),asy(i,2))];
end

K = size(eigenpairs,1);
xmesh = linspace(0,1,10000);
options = bvpset('RelTol',1e-8,'AbsTol',1e-10,'NMax',5000);

Y1 = cell(K,1);
Y2 = cell(K,1);

% Check solutions using bvp4c
for k = 1:K
    lambda1 = eigenpairs(k,1);
    lambda2 = eigenpairs(k,2);

    % Initial guess
    solinit = bvpinit(xmesh, [0; 1]);

    % Solve for y1
    sol1 = bvp4c(@(x,y) odefun_y1(x,y,lambda1,lambda2), ...
                 @bc_dirichlet, solinit, options);

    % Solve for y2
    sol2 = bvp4c(@(x,y) odefun_y2(x,y,lambda1,lambda2), ...
                 @bc_dirichlet, solinit, options);

    % Store eigenfunctions
    Y1{k} = deval(sol1, xmesh, 1);
    Y2{k} = deval(sol2, xmesh, 1);

    % Optional normalization
    Y1{k} = Y1{k} / norm(Y1{k});
    Y2{k} = Y2{k} / norm(Y2{k});
end

set(gca,'fontsize', 14);
subplot(2,1,1)
plot(xmesh,Y1{1},'LineWidth',.5)
subplot(2,1,2)
plot(xmesh,Y2{1},'LineWidth',.5)

function dydx = odefun_y1(x,y,lambda1,lambda2)
% y = [y1; y1']
dydx = zeros(2,1);
dydx(1) = y(2);
dydx(2) = -(lambda1*x + lambda2*(1-x))*y(1);
end

function dydx = odefun_y2(x,y,lambda1,lambda2)
% y = [y2; y2']
dydx = zeros(2,1);
dydx(1) = y(2);
dydx(2) = -(lambda1 + lambda2*x)*y(1);
end

function res = bc_dirichlet(ya,yb)
res = [ya(1); yb(1)];
end


function eg = testmpsl(x0,y0)
    % Set up boundary value problem
    % y_1''​(x)+{λ_1​x+λ_2​(1−x)}y_1​(x),    y_2''​(x)+{λ_1​+λ_2​x}*y_2​(x)​=0
    % [x,zt,A,B,C,G,kd,rd] = BDE2MEP(a,b,p,q,r,s,t,bc,N,opts) discretizes a two-parameter DE
    % p(x)y''(x) + q(x)y'(x) + r(x)y(x) = lambda s(x)y(x) + mu t(x)y(x)

    bc = [1 0; 1 0];
    n = 200;

    a1 = 0; 
    b1 = 1;
    p1 = 1;
    q1 = 0;
    r1 = @(x) x0*x + y0*(1-x);
    s1 = @(x) -x;
    t1 = @(x) x-1;
    [z1,A1,B1,C1,G1,k1,r1] = bde2mep(a1,b1,p1,q1,r1,s1,t1,bc,n);

    a2 = 0; 
    b2 = 1;
    p2 = 1;
    q2 = 0;
    r2 = @(y) x0 + y0*y;
    s2 = -1;
    t2 = @(y) -y;
    [z2,A2,B2,C2,G2,k2,r2] = bde2mep(a2,b2,p2,q2,r2,s2,t2,bc,n);


    % Set up opts
    opts.p = 50;
    opts.mom = 3;
    opts.block = 8;
    opts.h = 1;

    r = 50;

    P1 = zeros(size(A1,1),size(A1,2),2,2); P2 = P1;
    P1(:,:,1,1) = A1;
    P1(:,:,2,1) = B1*r;
    P1(:,:,1,2) = C1*r;
    P2(:,:,1,1) = A2;
    P2(:,:,2,1) = B2*r;
    P2(:,:,1,2) = C2*r;


    eg = biBeynLin(P1,P2,opts);
    if size(eg,1) > 0
        eg = eg.*r + [x0,y0];
    end

end