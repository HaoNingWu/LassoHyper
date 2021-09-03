% IMPORTANT! Please add the Chebfun Toolbox "chebfun-master" on to path 
% before running this demo

% Lasso hyperinterpolation on the interval

% Please note that our plotting are set for accommodating the plotting 
% requirments in our paper. If you choose to test another function, you may
% need to modify plotting setting.

clear, close all

% number of points and degree of polynomial
N = 300; L = 250; 


% mesh for plotting
xx = -1:0.001:1; 


% penalty parameter
mu = ones(1,L+1);


% basis index, 1 for Chebyshev polynomials of the 1st kind, 2 for Legendre
% polynomial
basis_idx = 2;
% generating corresponding nodes
switch basis_idx
    case 1
        [x,w] = jacpts(N+1,-.5,-.5); 
    case 2
        [x,w] = legpts(N+1);
end

% function to be approximated
example_idx = 4;
switch example_idx
    case 1
        GG = sin(pi*xx*5)./(pi*xx); GG(1001) = 5;
        G = sin(pi*x*5)./(pi*x);
        if mod(N+1,2) == 1
            G((N+2)/2) = 5;
        end
    case 2
        GG =  airy(40*xx); G= airy(40*x);    
    case 3    
        G = sin(25*x)+cos(16*x); GG = sin(25*xx)+cos(16*xx);
    case 4
        G = exp(-x.^2); GG = exp(-xx.^2);
end


% add noise
sigma = 0.15;
Y = G+normrnd(0,sigma,size(G));

% Generating matrix A
switch basis_idx
    case 1
        for l = 0:L
            for j = 0:N
                A(j+1,l+1) = cos(l*acos(x(j+1)))/sqrt(pi/2);
            end
        end
        A(:,1) = A(:,1)/sqrt(2);
    case 2
        for l = 0:L
            F = legpoly(l);
            A(:,l+1) = F(x)/sqrt(2/(2*l+1));
        end
end


% hyperinterpolation coefficients
beta0 = l1_beta(w,A,Y,0,mu); 

% L2 hyperinterpolation coefficients
LAMBDA = 10^(-1.5);
betatik = l2_beta(w,A,Y,LAMBDA,mu);


% filtered hyperinterpolation coefficients
for l = 0:1:L
    h(l+1) = Filter(L,l);
end
beta2 = h.*beta0;


% lasso hyperinterpolation coefficients
beta1 = l1_beta(w,A,Y,LAMBDA,mu);


% approximation polynomial on xx = -1:.001:1
p2 = zeros(2001,1); p1 = zeros(2001,1); p0 = zeros(2001,1);
for l = 0:L        
    switch basis_idx
        case 1
            if l == 0
                F = chebpoly(l,1)/sqrt(pi);
            else
                F = chebpoly(l,1)/sqrt(pi/2);
            end
        case 2
            F = legpoly(l)/sqrt(2/(2*l+1));
    end
    p2 = p2+beta2(l+1)*F(xx');
    p1 = p1+beta1(l+1)*F(xx');
    p0 = p0+betatik(l+1)*F(xx');
end
 
%% Figure 1
Color = parula(7);
fontsize_baseline = 20;
fontsize_baselinet = 30;
fontsize_baselinea = 15;
figure(1)
axes('position',[0.075 0.55 0.2 0.4]), 
plot(xx,GG,'linewidth',2,'color','k'), box on,...
    set(gca, 'fontsize', fontsize_baselinea),...
    xlabel('$x$','interpreter','latex', 'fontsize', fontsize_baseline),...%ylabel('$f(x)$','interpreter','latex', 'fontsize', fontsize_baseline),...
    title('exact function','interpreter','latex', 'fontsize', fontsize_baselinet),...
     grid on,...
     set(gca, 'XMinorGrid', 'off'), set(gca, 'YMinorGrid', 'off'),axis([-1,1,-.4,1.6]),
axes('position',[0.075 0.05 0.2 0.4]), 
plot(x,Y,'linewidth',2,'color','k'),box on,...
    set(gca, 'fontsize', fontsize_baselinea),...
    xlabel('$x$','interpreter','latex', 'fontsize', fontsize_baseline),...%ylabel('$f(x)$','interpreter','latex', 'fontsize', fontsize_baseline),...
    title('noisy function','interpreter','latex', 'fontsize', fontsize_baselinet),...
     grid on,...
     set(gca, 'XMinorGrid', 'off'), set(gca, 'YMinorGrid', 'off'),axis([-1,1,-.4,1.6]),
axes('position',[0.3 0.55 0.2 0.4]),
plot(xx,p0,'linewidth',2,'color','k'), box on,...
    set(gca, 'fontsize', fontsize_baselinea),...
    xlabel('$x$','interpreter','latex', 'fontsize', fontsize_baseline),...%ylabel('$f(x)$','interpreter','latex', 'fontsize', fontsize_baseline),...
    title('$\mathcal{T}_Lf$','interpreter','latex', 'fontsize', fontsize_baselinet),  grid on,...
     set(gca, 'XMinorGrid', 'off'), set(gca, 'YMinorGrid', 'off'),axis([-1,1,-.4,1.6]),
axes('position',[0.3 0.05 0.2 0.4]), 
plot(xx,abs(GG-p0'),'linewidth',2,'color','k'),set(gca, 'fontsize', fontsize_baselinea),...
    xlabel('$x$','interpreter','latex', 'fontsize', fontsize_baseline),...%ylabel('Absolute error', 'interpreter','latex','fontsize', fontsize_baseline),...
    title('error','interpreter','latex','fontsize', fontsize_baselinet),box on, grid on,...
     set(gca, 'XMinorGrid', 'off'), set(gca, 'YMinorGrid', 'off'),axis([-1,1,0,1.5]),
axes('position',[0.525 0.55 0.2 0.4]),
plot(xx,p2,'linewidth',2,'color','k'),box on,set(gca, 'fontsize', fontsize_baselinea), ...
     xlabel('$x$','interpreter','latex', 'fontsize', fontsize_baseline),...%ylabel('$f(x)$','interpreter','latex', 'fontsize', fontsize_baseline),...
     title('$\mathcal{F}_Lf$','interpreter','latex', 'fontsize', fontsize_baselinet),...
     grid on,...
     set(gca, 'XMinorGrid', 'off'), set(gca, 'YMinorGrid', 'off'),axis([-1,1,-.4,1.6]),
axes('position',[0.525 0.05 0.2 0.4]), set(gca, 'fontsize', fontsize_baselinea), 
plot(xx,abs(GG-p2'),'linewidth',2,'color','k'),set(gca, 'fontsize', fontsize_baseline),...
    xlabel('$x$','interpreter','latex', 'fontsize', fontsize_baseline),...%ylabel('Absolute error','interpreter','latex', 'fontsize', fontsize_baseline),...
    title('error','interpreter','latex', 'fontsize', fontsize_baselinet),box on,axis([-1,1,0,1.5]), grid on,...
     set(gca, 'XMinorGrid', 'off'), set(gca, 'YMinorGrid', 'off')
axes('position',[0.75 0.55 0.2 0.40]),
plot(xx,p1,'linewidth',2,'color','k'),box on,set(gca, 'fontsize', fontsize_baselinea), ...
     xlabel('$x$','interpreter','latex', 'fontsize', fontsize_baseline),...%ylabel('$f(x)$','interpreter','latex', 'fontsize', fontsize_baseline),...
     title('$\mathcal{L}_L^{\lambda}f$ with $\lambda=10^{-1}$','interpreter','latex', 'fontsize', fontsize_baselinet),...
     grid on,...
     set(gca, 'XMinorGrid', 'off'), set(gca, 'YMinorGrid', 'off'),axis([-1,1,-.4,1.6]),
axes('position',[0.75 0.05 0.2 0.40]), 
plot(xx,abs(GG-p1'),'linewidth',2,'color','k'),set(gca, 'fontsize', fontsize_baselinea), xlabel('$x$','interpreter','latex', 'fontsize', fontsize_baseline),...%ylabel('Absolute error','interpreter','latex', 'fontsize', fontsize_baseline),...
    title('error','interpreter','latex', 'fontsize', fontsize_baselinet),box on,axis([-1,1,0,1.5]),grid on,...
     set(gca, 'XMinorGrid', 'off'), set(gca, 'YMinorGrid', 'off')

 
 
 