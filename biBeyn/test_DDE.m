%%%Test of biBeyn on analytic problem 

% Set up problem
sols2 = [];
sols3 = [];
h2 = 3:-.1:1;
n = 2;

rng(2);
A1 = randn(n);
A2 = randn(n);
A0 = randn(n);

xest = 2.1;
yest = 2.7;

for i=1:numel(h2)

a = h2(i);

P1 = @(x,y) -eye(n).*1i*x + A0 + A1.*exp(-1i.*x.*y) + A2.*exp(-1i.*x.*h2(i));
P2 = @(x,y) eye(n).*1i.*x + A0 + A1.*exp(1i*x*y) + A2.*exp(1i*x*h2(i));

tmp = critdelays(A0,A1,A2,a,xest,yest);
if numel(tmp) > 0
        sols2 = [sols2; tmp(1,2), a];
end

sols = tmp;
xest = tmp(1,1);
yest = tmp(1,2);

% local check
F = @(z) [det(P1(z(1),z(2))) det(P2(z(1),z(2)))];
solloc = fsolve(F,tmp);
sols3 = [sols3; solloc(2), a];
end

% plot
figure;
set(gca,'fontsize', 14);
plot(sols3(:,1), sols3(:,2), 'bo-', 'LineWidth', 1.5, 'MarkerSize', 6);
grid on;



function eg = critdelays(A0,A1,A2,h2,xest,hest)
    % finds critical delays for a given set of matrices
    % characteristic equation and conjugate
    rng(3);
    n = size(A0,1);
    P1 = @(x,y) -eye(n).*1i*x + A0 + A1.*exp(-1i.*x.*y) + A2.*exp(-1i.*x.*h2);
    P2 = @(x,y) eye(n).*1i.*x + A0 + A1.*exp(1i*x*y) + A2.*exp(1i*x*h2);

   dP1dx = @(x,y) ...
    -1i*eye(n) ...
  + A1.*(-1i*y).*exp(-1i*x*y) ...
  + A2.*(-1i*h2).*exp(-1i*x*h2);


    dP1dy = @(x,y) ...
    A1.*(-1i*x).*exp(-1i*x*y);


    dP2dx = @(x,y) ...
    1i*eye(n) ...
  + A1.*(1i*y).*exp(1i*x*y) ...
  + A2.*(1i*h2).*exp(1i*x*h2);


    dP2dy = @(x,y) ...
    A1.*(1i*x).*exp(1i*x*y);

    % settings for biBeyn
    opts.block = 3;
    opts.mom = 20;
    opts.p = 200;

    % r = 2;
    eg = [];
    x = xest; y = hest;
    zx = 1; zy = 1;
    sols = biBeyn(@(z1,z2) P1(x+z1/zx,y+z2/zy),@(z1,z2) P2(x+z1/zx,y+z2/zy),...
        @(z1,z2) dP1dx(x+z1/zx,y+z2/zy),@(z1,z2) dP1dy(x+z1/zx,y+z2/zy),...
        @(z1,z2) dP2dx(x+z1/zx,y+z2/zy),@(z1,z2) dP2dy(x+z1/zx,y+z2/zy),opts);
    if numel(sols) > 0
        sols(:,1) = sols(:,1)./zx+x; sols(:,2) = sols(:,2)./zy+y;
        eg = sols;
    end
end