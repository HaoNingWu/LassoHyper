% IMPORTANT! Please add the Chebfun Toolbox "sphere_approx_toolbox_v3.0" 
% onto path before running this demo

% Lasso hyperinterpolation on the interval

% Please note that our plotting are set for accommodating the plotting 
% requirments in our paper. If you choose to test another function, you may
% need to modify plotting setting.

clear ; close all

L = 25;

    % function information
    funtxt_ce = {'nrWend_k0','nrWend_k1','nrWend_k2','nrWend_k3','nrWend_k4'};
    i_f = 3;
    funtxt = funtxt_ce{i_f};
    switch funtxt
        case 'nrWend_k0'
            rbf_k = 0;
        case 'nrWend_k1'
            rbf_k = 1;
        case 'nrWend_k2'
            rbf_k = 2;
        case 'nrWend_k3'
            rbf_k = 3;
        case 'nrWend_k4'
            rbf_k = 4;
    end
    switch funtxt
         case {'nrWend_k0','nrWend_k1','nrWend_k2','nrWend_k3','nrWend_k4'}
             func = @rbf_nr;
    end
    k = rbf_k;
    delta = (3*k+3)*gamma(k+1/2)/(2*gamma(k+1));



    model_parameter.func = func;



    % Validation point set
    Xt = get_Xt( );
    ft = func(Xt',rbf_k);
    [ Yt ] = get_Yt( L, Xt );
    model_parameter.func_delta = 0;


    LAMBDA = [0 10^(-3.5) 10^(-2.5) 10^(-2)];



    noise = 0.015; % Gaussian noise
    r = 0.02; % impulse noise


    t_now = 2*L;
    % degree of point set and polynomial
    model_parameter.t = t_now;
    model_parameter.L = L;

    X_k = loadStd( model_parameter.t, (model_parameter.t+1)^2 );
    
    % generating function
    f = func(X_k,rbf_k);
    
    % add noise
    f = f+ normrnd(0,noise,(model_parameter.t+1)^2,1)+r*binornd(1,r,(model_parameter.t+1)^2,1);
    
    % hyperinterpolation coefficients
    model_parameter.beta = 1;  
    model_parameter.reg_type = 3;
    model_parameter.reg_parameter.s = 1;
    model_parameter.lambda = 0;
    res_now_0 = Solver( 'L2L2', X_k, f, model_parameter);
    alpha0 =  res_now_0.alpha;
    
    % Tikhonov hyperinterpolation coefficients with Laplace-Beltrami
    % operator
    model_parameter.beta = 1;   
    model_parameter.reg_type = 3;
    model_parameter.reg_parameter.s = 1;
    model_parameter.lambda = LAMBDA(2);
    res_now_1 = Solver( 'L2L2', X_k, f, model_parameter);
    alpha1 =  res_now_1.alpha;    
    
    % Lasso hyperinterpolation coefficients
    for i =1:length(alpha0)
        alpha2(i,1) =  wthresh(alpha0(i),'s',LAMBDA(3));
    end
    
    % filtered hyperinterpolation coefficients
    model_parameter.reg_parameter.type_filter = 2; 
    model_parameter.reg_type = 2;     
    model_parameter.lambda = 1;
    res_now_4 = Solver( 'L2L2', X_k, f, model_parameter);
    alpha4 =  res_now_4.alpha;





    % approximation polynomials
    L = sqrt(length(alpha2))-1;
    Yt_eqp = get_Yt( L, Xt );
    pt_eqptik = Yt_eqp * alpha1;
    pt_eqplasso = Yt_eqp * alpha2;
    pt_eqpfiltered = Yt_eqp * alpha4;




    % scale the plotting for better visual comparison 
    Xt = Xt';
    x = Xt(:,1); y = Xt(:,2); z = Xt(:,3);
    tri = convhull([x y z]);


    C1 = ft;
    C3 = pt_eqptik;
    C4 = pt_eqplasso;
    C5 = pt_eqpfiltered;


    [Fmax, imax] = max(C1);
    [Fmin, imin] = min(C1);
    scale = 0.5;
    FS = 1 + (scale/(Fmax-Fmin))*(C1-Fmin);

    fontsize_baselinet = 30;
    marksize = 80;
    color = lines(4);
    count_hx = 2;

    % plotting

    axes('position',[0,0.55,0.235,0.48]) 
    fg = trisurf(tri,x.*FS, y.*FS, z.*FS, C1,'facecolor','interp');
    set(fg,'EdgeColor', 'none'); 
    title('exact function','interpreter','latex','fontsize',fontsize_baselinet,'position',[-.8,.9,1])
    colormap(jet(255));
    view(39,46), axis vis3d, axis equal tight, colorbar('south'), caxis([1.42,1.57])
    axis off
   
    
    axes('position',[0.24,0.55,0.235,0.48])
    [Fmax, imax] = max(C3);
    [Fmin, imin] = min(C3);
    scale = 0.5;
    FS = 1 + (scale/(Fmax-Fmin))*(C3-Fmin);
    fg = trisurf(tri,x.*FS, y.*FS, z.*FS, C3,'facecolor','interp');
    set(fg,'EdgeColor', 'none'); 
    title('$\mathcal{T}_Lf$ with $\lambda=10^{-3.5}$','interpreter','latex','fontsize',fontsize_baselinet,'position',[-.8,.9,1])
    colormap(jet(255));
    view(39,46), axis vis3d, axis equal tight, colorbar('south'),caxis([1.42,1.57])
    axis off 
   
    axes('position',[0.72,0.55,0.235,0.48])
    [Fmax, imax] = max(C4);
    [Fmin, imin] = min(C4);
    scale = 0.5;
    FS = 1 + (scale/(Fmax-Fmin))*(C4-Fmin);
    fg = trisurf(tri,x.*FS, y.*FS, z.*FS, C4,'facecolor','interp');
    set(fg,'EdgeColor', 'none'); 
    title('$\mathcal{L}_L^{\lambda}f$ with $\lambda=10^{-2.5}$','interpreter','latex','fontsize',fontsize_baselinet,'position',[-.8,0.9,1])
    colormap(jet(255));
    view(39,46), axis vis3d, axis equal tight, colorbar('south'),caxis([1.42,1.57])
    axis off  
 
   
    axes('position',[0.48,0.55,0.235,0.48])
    [Fmax, imax] = max(C5);
    [Fmin, imin] = min(C5);
    scale = 0.5;
    FS = 1 + (scale/(Fmax-Fmin))*(C5-Fmin);
    fg = trisurf(tri,x.*FS, y.*FS, z.*FS, C5,'facecolor','interp');
    set(fg,'EdgeColor', 'none'); 
    title('$\mathcal{F}_Lf$','interpreter','latex','fontsize',fontsize_baselinet,'position',[-.8,.9,1])
    colormap(jet(255));
    view(39,46), axis vis3d, axis equal tight, colorbar('south'),caxis([1.42,1.57])
    axis off
   

     axes('position',[0.24,0.05,0.235,0.48])
    [Fmax, imax] = max(abs(C3-C1));
    [Fmin, imin] = min(abs(C3-C1));
    scale = 0.5;
    FS = 1 + (scale/(Fmax-Fmin))*(abs(C3-C1)-Fmin);
    fg = trisurf(tri,x.*FS, y.*FS, z.*FS, abs(C3-C1),'facecolor','interp');
    set(fg,'EdgeColor', 'none'); 
    title('error','interpreter','latex','fontsize',fontsize_baselinet,'position',[-.8,.9,1])
    colormap(jet(255));
    view(39,46), axis vis3d, axis equal tight, colorbar('south'),caxis([0,0.016])
    axis off  

 
    axes('position',[0.48,0.05,0.235,0.48])
    [Fmax, imax] = max(abs(C5-C1));
    [Fmin, imin] = min(abs(C5-C1));
    scale = 0.5;
    FS = 1 + (scale/(Fmax-Fmin))*(abs(C5-C1)-Fmin);
    fg = trisurf(tri,x.*FS, y.*FS, z.*FS, abs(C5-C1),'facecolor','interp');
    set(fg,'EdgeColor', 'none'); 
    title('error','interpreter','latex','fontsize',fontsize_baselinet,'position',[-.8,.9,1])
    colormap(jet(255));
    view(39,46), axis vis3d, axis equal tight, colorbar('south'),caxis([0,0.016])
    axis off
   
   
    axes('position',[0.72,0.05,0.235,0.48])
    [Fmax, imax] = max(abs(C4-C1));
    [Fmin, imin] = min(abs(C4-C1));
    scale = 0.5;
    FS = 1 + (scale/(Fmax-Fmin))*(abs(C4-C1)-Fmin);
    fg = trisurf(tri,x.*FS, y.*FS, z.*FS, abs(C4-C1),'facecolor','interp');
    set(fg,'EdgeColor', 'none'); 
    title('error','interpreter','latex','fontsize',fontsize_baselinet,'position',[-.8,.9,1])
    colormap(jet(255));
    view(39,46), axis vis3d, axis equal tight, colorbar('south'),caxis([0,0.016])
    axis off
   
    
    
    % plot noisy function
    Yt_eqp = get_Yt( L, X_k );
    Xt = X_k;
    x = Xt(:,1); y = Xt(:,2); z = Xt(:,3);
    tri = convhull([x y z]);
    C2 = f;
    
    axes('position',[0,0.05,0.235,0.48])
    [Fmax, imax] = max(C2);
    [Fmin, imin] = min(C2);
    scale = 0.5;
    FS = 1 + (scale/(Fmax-Fmin))*(C2-Fmin);
    fg = trisurf(tri,x.*FS, y.*FS, z.*FS, C2,'facecolor','interp');
    set(fg,'EdgeColor', 'none'); 
    title('noisy function','interpreter','latex','fontsize',fontsize_baselinet,'position',[-.8,.9,1])
    colormap(jet(255));
    view(39,46), axis vis3d, axis equal tight, colorbar('south'),caxis([1.42,1.57]),
    axis off

