% test_ARMA11_2  demo for ARMA11 using multiBeyn
%
% We find critial points for the ARMA(1,1) model
%
% Based on
% M.E.Hochstenbach, T.Kosir, B.Plestenjak: On the 
% solution of rectangular multiparameter eigenvalue problems, arXiv 2212.01867

% MultiParEig toolbox
% B. Plestenjak, University of Ljubljana
% FreeBSD License, see LICENSE.txt

% Contour solver from
% multiBeyn
% E. Graf and A. Townsend, Cornell University
% 2026

rng(3);
%% 2nd test - match multipareig on a short sequence

y = [2.4130 1.0033 1.2378 -0.72191 -0.81745 -2.2918 0.18213 0.073557 0.55248 2.0180 2.6593 1.1791];
y = y(:);

opts = [];
opts.rtol = 1e-3;
opts.xbd = [-1 1];
opts.ybd = [-1 1];
[points,val,err,cand] = arma11(y,opts);
[points2,val2,err2,cand2] = ARMA11(y,opts);


N = length(y);
M = 800;
aset = linspace(-2,2,M);
gset = linspace(-2,2,M);
Z = zeros(M,M);
parfor i = 1:M
    for j = 1:M
        TC = diag(gset(j)*ones(N-1,1))+diag(ones(N-2,1),1); TC(N-1,N) = 1;
        TA = diag(aset(i)*ones(N-1,1))+diag(ones(N-2,1),1); TA(N-1,N) = 1;
        e = pinv(TC)*TA*y;
        Z(j,i) = norm(e)^2;
    end
end

contour(aset,gset,log(Z),60,'LineWidth',3);
hold on
alpha = real(points(:,1));
gamma = real(points(:,2));
alpha2 = real(points2(:,1));
gamma2 = real(points2(:,2));
plot(alpha,gamma,'bo','MarkerSize',12,'Markerfacecolor','b','LineWidth',1.5);
plot(alpha2,gamma2,'ro','MarkerSize',8,'Markerfacecolor','r','LineWidth',1.5);
hold off
xlabel('\alpha_1','FontSize',14)
ylabel('\gamma_1','FontSize',14)