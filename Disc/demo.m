% Lasso hyperinterpolation on the disc

% Please note that our plotting are set for accommodating the plotting 
% requirments in our paper. If you choose to test another function, you may
% need to modify plotting setting.

clear all, close all

% Input problem parameters
degree = 16; 
num_quad = 136; 
index_f = 5;

% Initialization
nd = degree; 
nq = num_quad;
order =((nd+1)*(nd+2))/2;

    
lambda = 10^(-1.5);
noise = 3.5;


% Computing coefficients
% alpha: hyperinterpolation coefficients
% alphafiltered: filtered hyperinterpolation coefficients
[alpha,alphafiltered,noi] = hyperdisc(nd,index_f,nq,noise);

% Lasso hyperinterpolation coefficients
for i = 1:order
    if alpha(i)>lambda
        beta(i) = alpha(i)-lambda;
    elseif alpha(i)<-lambda
        beta(i) = alpha(i)+lambda;
    else
        beta(i) = 0;
    end
end

for i = 1:order
%     beta(i) = wthresh(alpha(i),'h',lambda)
    tik(i) = alpha(i)/(1+lambda*(1/Filter(order,i-1))^2);
end

% generating approximation polynomials
pr = 36;
r = 0:1/pr:1;
th = 0:pi/pr:2*pi;
[R,TH] = meshgrid(r,th);
X = R.*cos(TH); Y = R.*sin(TH);
tik_sol = eval_poly3(X,Y,nd,tik);
filter_sol = eval_poly3(X,Y,nd,alphafiltered);
lasso_sol = eval_poly3(X,Y,nd,beta);



true_sol = true_solution(X,Y,index_f);


% Plotting
fontsize_baselinea = 15;
fontsize_baseline = 20;
fontsize_baselinet = 30;

axes('position',[0.075 0.55 0.2 0.4]), 
mesh(X,Y,true_sol,'edgecolor','k'), set(gca, 'fontsize', fontsize_baselinea), box on,...
    xlabel('$x$','interpreter','latex', 'fontsize', fontsize_baseline),ylabel('$y$','interpreter','latex', 'fontsize', fontsize_baseline),...%ylabel('$f(x)$','interpreter','latex', 'fontsize', fontsize_baseline),...
    title('exact function','interpreter','latex', 'fontsize', fontsize_baselinet),...
      grid on, axis([-1,1,-1,1,-.5,1.5]),

axes('position',[0.3 0.55 0.2 0.4]),
mesh(X,Y,tik_sol,'edgecolor','k'),set(gca, 'fontsize', fontsize_baselinea),box on,...
    xlabel('$x$','interpreter','latex', 'fontsize', fontsize_baseline),ylabel('$y$','interpreter','latex', 'fontsize', fontsize_baseline),...%ylabel('$f(x)$','interpreter','latex', 'fontsize', fontsize_baseline),...
    title('$\mathcal{T}_Lf$','interpreter','latex', 'fontsize', fontsize_baselinet),  grid on,axis([-1,1,-1,1,-.5,1.5]),

axes('position',[0.3 0.05 0.2 0.4]), 
error = abs(true_sol - tik_sol);
mesh(X,Y,error,'edgecolor','k'),set(gca, 'fontsize', fontsize_baselinea), xlabel('$x$','interpreter','latex', 'fontsize', fontsize_baseline), ylabel('$y$','interpreter','latex', 'fontsize', fontsize_baseline),...%ylabel('Absolute error', 'interpreter','latex','fontsize', fontsize_baseline),...
    title('error','interpreter','latex','fontsize', fontsize_baselinet),box on,  grid on,...
%      set(gca, 'XMinorGrid', 'off'), set(gca, 'YMinorGrid', 'off'),axis([-1,1,0,1.5]),
 axis([-1,1,-1,1,0,0.5])
 
 
axes('position',[0.525 0.55 0.2 0.4]),
mesh(X,Y,filter_sol,'edgecolor','k'), set(gca, 'fontsize', fontsize_baselinea),box on,...
     xlabel('$x$','interpreter','latex', 'fontsize', fontsize_baseline),ylabel('$y$','interpreter','latex', 'fontsize', fontsize_baseline),...%ylabel('$f(x)$','interpreter','latex', 'fontsize', fontsize_baseline),...
     title('$\mathcal{F}_Lf$','interpreter','latex', 'fontsize', fontsize_baselinet),...
     grid on,...
%      set(gca, 'XMinorGrid', 'off'), set(gca, 'YMinorGrid', 'off'),
axis([-1,1,-1,1,-.5,1.5]),

axes('position',[0.525 0.05 0.2 0.4]), 
error = abs(true_sol - filter_sol);
mesh(X,Y,error,'edgecolor','k'),set(gca, 'fontsize', fontsize_baselinea),xlabel('$x$','interpreter','latex', 'fontsize', fontsize_baseline),ylabel('$y$','interpreter','latex', 'fontsize', fontsize_baseline),...%ylabel('Absolute error','interpreter','latex', 'fontsize', fontsize_baseline),...
    title('error','interpreter','latex', 'fontsize', fontsize_baselinet),box on,  grid on,...
     set(gca, 'XMinorGrid', 'off'), set(gca, 'YMinorGrid', 'off'), axis([-1,1,-1,1,0,0.5])

 axes('position',[0.75 0.55 0.2 0.4]),
mesh(X,Y,lasso_sol,'edgecolor','k'),set(gca, 'fontsize', fontsize_baselinea),box on,...
     xlabel('$x$','interpreter','latex', 'fontsize', fontsize_baseline),ylabel('$y$','interpreter','latex', 'fontsize', fontsize_baseline),...%ylabel('$f(x)$','interpreter','latex', 'fontsize', fontsize_baseline),...
     title('$\mathcal{L}_L^{\lambda}f$ with $\lambda=10^{-1.5}$','interpreter','latex', 'fontsize', fontsize_baselinet),...
      grid on,...
%      set(gca, 'XMinorGrid', 'off'), set(gca, 'YMinorGrid', 'off'),
axis([-1,1,-1,1,-.5,1.5]),

axes('position',[0.75 0.05 0.2 0.4]), 
error = abs(true_sol - lasso_sol);
mesh(X,Y,error,'edgecolor','k'),set(gca, 'fontsize', fontsize_baselinea),xlabel('$x$','interpreter','latex', 'fontsize', fontsize_baseline),ylabel('$y$','interpreter','latex', 'fontsize', fontsize_baseline),...%ylabel('Absolute error','interpreter','latex', 'fontsize', fontsize_baseline),...
    title('error','interpreter','latex', 'fontsize', fontsize_baselinet),box on, grid on,...
     set(gca, 'XMinorGrid', 'off'), set(gca, 'YMinorGrid', 'off'),  axis([-1,1,-1,1,0,0.5])





% When the number of radial points pr = 36, the noisy plot is a mess! We
% reduce the number to 18 but with the same level of noise, trying to
% display the perturbed function
pr = 18;
r = 0:1/pr:1;
th = 0:pi/pr:2*pi;
[R,TH] = meshgrid(r,th);
X = R.*cos(TH); Y = R.*sin(TH);
true_sol = true_solution(X,Y,index_f);
 
axes('position',[0.075 0.05 0.2 0.4]), 
mesh(X,Y,true_sol+noise*randn(size(X)),'edgecolor','k'),set(gca, 'fontsize', fontsize_baselinea),box on,...
    xlabel('$x$','interpreter','latex', 'fontsize', fontsize_baseline),ylabel('$y$','interpreter','latex', 'fontsize', fontsize_baseline),...%ylabel('$f(x)$','interpreter','latex', 'fontsize', fontsize_baseline),...
    title('noisy function','interpreter','latex', 'fontsize', fontsize_baselinet),...
      grid on,...
axis([-1,1,-1,1,-10,10]),






function sum = eval_poly3(x,y,n,alpha)

% Evaluate the polynomial of degree n that is expressed as a sum of
% orthonormal ridge polynomials [each of which is also multiplied by the
% factor of (1-x^2-y^2)], with the coefficients given in alpha.

order = (n+1)*(n+2)/2;
sum = zeros(size(x));
for m=1:order
    [km,jm] = compute_KJ(m);
    sum = sum + alpha(m)*ridge_poly(x,y,km,jm);
end

end 


%====================================================================

function [K,J] = compute_KJ(n)

% Compute the degree K and the index J for the ridge
% polynomial of lexicographic index n.

K = floor((-1+sqrt(1+8*n))/2-10*eps);
J = n - (K*(K+1))/2 - 1;

end 


% ====================================================================

function ans=ridge_poly(x,y,n,k)

% This evaluates the ridge orthonormal polynomial of
% degree n and index k.

arg = k*pi/(n+1);
c_x = cos(arg); c_y = sin(arg);
z = c_x*x + c_y*y;
ans = (1/sqrt(pi))*cheb2_recurs(z,n);

end

% =========================================================

function U=cheb2_recurs(x,n)

% This evaluates at x the Chebyshev polynomial of the second
% kind of degree n using the triple recursion relation for those
% polynomials.

if n==0
    U = ones(size(x));
    return
end

if n==1
    U = 2*x;
    return
end

% n > 1. Use recursion to define U.

[L1,L2]=size(x);
S = zeros(L1,L2,n+1);       % Dimension

S(:,:,1) = ones(size(x));
x2 = 2*x;
S(:,:,2) = x2;
for k=2:n
    S(:,:,k+1) = x2.*S(:,:,k) - S(:,:,k-1);
end

U = S(:,:,n+1);

end 

% ====================================================================

function ans = true_solution(x,y,index_f)
switch index_f
case 1
    ans = 10*(1-(x.^2+y.^2));
case 2
    ans = 3 + (x.^2+y.^2);
case 3
    ans = 4 - (1 - (x.^2+y.^2)).^2;
case 5
    ans = (1-(x.^2+y.^2)).*exp(x.*cos(y));
case 6
    ans = 1./(1+(x.^2+y.^2));
end

end % true_solution



