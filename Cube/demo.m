
% ----------------------------------------
%  This M-file allows to compute the 
%  interpolant, either using non-tensor
% Lasso hyperinterpolation in the cube

% Please note that our plotting are set for accommodating the plotting 
% requirments in our paper. If you choose to test another function, you may
% need to modify plotting setting.

clear, close all

n = 50; % degree

f_idx = 5;

sigma = 0.2; % level of noise
lambda =10^(-2.5);



        
         
[err_I,errest_I,ERRORTIKtmp,LV_I,LCalfa_I,Lcputime_I,Calfa,fname,L,Fhyper] = CX_Interpolation(n,f_idx,lambda,sigma);
[err_Ilasso, errest_Ilasso,ERRORLASSOtmp,LV_lasso, LCalfa_Ilasso, Lcputime_Ilasso,Calfalasso,fname,F] = CX_Interpolationlasso(n,f_idx,lambda,sigma);
[err_Ifilter, errest_Ifilter,ERRORFILtmp,LV_filter, LCalfa_Ifilter, Lcputime_Ifilter,Calfafilter,fname1,Ffiltered] = CX_Interpolationfiltered(n,f_idx,sigma);








dis = 0.1;
[xx,yy,zz] = meshgrid(-1:dis:1);
nnn = 2/dis+1;

f_targets = @(x,y,z) exp(-1./(x.^2+y.^2+z.^2));
Vexact = f_targets(xx,yy,zz);
Vnoisy = Vexact;
Vnoisy(find(Vnoisy~=0)) = noisegen(Vnoisy(find(Vnoisy~=0)),10);
Vnoisy = reshape(Vnoisy,[nnn,nnn,nnn]);


Vhyper =Fhyper(xx,yy,zz);

Vlasso =F(xx,yy,zz);

Vfiltered =Ffiltered(xx,yy,zz);



% plotting

fontsize_baselinea = 15;
fontsize_baseline = 20;
fontsize_baselinet = 30;
marksize = 80;

colormap(jet)


xslice = [-.25,.5,1]; 
yslice = [0,1]; 
zslice = [-1,0];


axes('position',[0.02,0.55,0.2,0.4]) 
slice(xx,yy,zz,Vexact,xslice,yslice,zslice),set(gca, 'fontsize', fontsize_baselinea),...
title('exact function','interpreter','latex', 'fontsize', fontsize_baselinet),...
grid on, set(gca, 'XMinorGrid', 'off'), set(gca, 'YMinorGrid', 'off'),  view(-36,15), colorbar('eastoutside')

axes('position',[0.02,0.05,0.2,0.4])
slice(xx,yy,zz,Vnoisy,xslice,yslice,zslice),set(gca, 'fontsize', fontsize_baselinea),...
title('noisy function','interpreter','latex', 'fontsize', fontsize_baselinet),...
grid on, set(gca, 'XMinorGrid', 'off'), set(gca, 'YMinorGrid', 'off'),  view(-36,15), colorbar('eastoutside'),caxis([0,0.7])

axes('position',[0.26,0.55,0.2,0.4])
slice(xx,yy,zz,Vhyper,xslice,yslice,zslice),set(gca, 'fontsize', fontsize_baselinea),...
title('$\mathcal{T}_Lf$','interpreter','latex', 'fontsize', fontsize_baselinet),...
grid on, set(gca, 'XMinorGrid', 'off'), set(gca, 'YMinorGrid', 'off'), view(-36,15), colorbar('eastoutside'), caxis([0,.7])


axes('position',[0.26,0.05,0.2,0.4])
slice(xx,yy,zz,abs(Vhyper-Vexact),xslice,yslice,zslice),set(gca, 'fontsize', fontsize_baselinea),...
title('error','interpreter','latex', 'fontsize', fontsize_baselinet),...
grid on, set(gca, 'XMinorGrid', 'off'), set(gca, 'YMinorGrid', 'off'),  view(-36,15), colorbar('eastoutside'), caxis([0,0.25])




axes('position',[0.5,0.55,0.2,0.4])
slice(xx,yy,zz,Vfiltered,xslice,yslice,zslice),set(gca, 'fontsize', fontsize_baselinea),...
title('$\mathcal{F}_Lf$','interpreter','latex', 'fontsize', fontsize_baselinet),...
grid on, set(gca, 'XMinorGrid', 'off'), set(gca, 'YMinorGrid', 'off'), view(-36,15),colorbar('eastoutside'),caxis([0,.7])
 
axes('position',[0.5,0.05,0.2,0.4])
slice(xx,yy,zz,abs(Vfiltered-Vexact),xslice,yslice,zslice),set(gca, 'fontsize', fontsize_baselinea),... 
title('error','interpreter','latex', 'fontsize', fontsize_baselinet),...
grid on, set(gca, 'XMinorGrid', 'off'), set(gca, 'YMinorGrid', 'off'),  view(-36,15),colorbar('eastoutside'), caxis([0,0.25])


axes('position',[0.74,0.55,0.2,0.4])
slice(xx,yy,zz,Vlasso,xslice,yslice,zslice),set(gca, 'fontsize', fontsize_baselinea),...
title('$\mathcal{L}_L^{\lambda}f$ with $\lambda=10^{-2.5}$','interpreter','latex', 'fontsize', fontsize_baselinet),...
grid on,set(gca, 'XMinorGrid', 'off'), set(gca, 'YMinorGrid', 'off'),  view(-36,15), colorbar('eastoutside'),caxis([0,.7])




axes('position',[0.74,0.05,0.2,0.4])
slice(xx,yy,zz,abs(Vlasso-Vexact),xslice,yslice,zslice),set(gca, 'fontsize', fontsize_baselinea),... 
title('error','interpreter','latex', 'fontsize', fontsize_baselinet),...
grid on, set(gca, 'XMinorGrid', 'off'), set(gca, 'YMinorGrid', 'off'),  view(-36,15),colorbar('eastoutside'),caxis([0,0.25])







