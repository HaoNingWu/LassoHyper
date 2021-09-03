function [alpha,alphafiltered,noi] = ...
              hyperdisc(degree,index_f,num_quad,noise)

% Computing hyperinterpolation and filtered hyperinterpolation coefficients

% Input parameters:
%   degree:     The degree of the polynomial approximation being sought.
%   index_f:    The index of the function f(x,y,u) defining the Poisson
%               nonlinear equation.
%   num_quad:   The index for the Gaussian quadrature being used.
%   error_tol:  The error tolerance for the Newton's method iteration.
              
% Initialization
nd = degree;
nq = num_quad;
order = ((nd+1)*(nd+2))/2;


% Initialize points and weights for quadrature over the disk.
hq = 2*pi/(2*nq+1);
t = hq*[1:2*nq+1];
cos_t = cos(t);
sin_t = sin(t);
[r_nodes,r_wgts] = gl_nodes_wts(0,1,nq+1);
r_wgts = r_nodes.*r_wgts;
[R,Th] = meshgrid(r_nodes,t);
X = R.*cos(Th); Y = R.*sin(Th);
clear R; clear Th;

h_vec = hq*ones(1,2*nq+1);
[Wradial,Wtheta] = meshgrid(r_wgts,h_vec);
weight = Wradial.*Wtheta;
disp(['area = ',num2str(sum(sum(weight)))])     %Print integral check
clear Wradial; clear Wtheta;


% data for approximation
fv = f_given(X,Y,index_f);
noi = noise*binornd(1,0.5,size(X)).*(1-2*rand(size(X)));
fv = fv+noi;

for n=1:order
    [km,jm] = compute_KJ(n); 
    phi_m = ridge_poly(X,Y,km,jm);
    alpha(n) = sum(sum(weight.*phi_m.*fv));
    alphafiltered(n)=Filter(nd,km)*sum(sum(weight.*phi_m.*fv));
end
    
    
end
    
%====================================================================

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
sum = sum.*(1 - (x.^2 + y.^2));

end % eval_poly3

%====================================================================

function [K,J] = compute_KJ(n)

% Compute the degree K and the index J for the ridge
% polynomial of lexicographic index n.

K = floor((-1+sqrt(1+8*n))/2-10*eps);
J = n - (K*(K+1))/2 - 1;

end % compute_KJ

%====================================================================

function value = Gram2(d1,k1,d2,k2)

% This performs numerical integration over the unit disk of the 
% products of the special polynomials Psi(x,y;d1,k1) 
% and Psi(x,y;d2,k2). The polynomial Psi(x,y;d,k) is of degree d 
% and index k, and it is defined by
%   Psi(x,y;d,k) = Laplacian[(1-x^2-y^2)Phi(x,y;d,k)]
% where Phi(x,y;d,k) is the ridge polynomial of degree d and
% index k.

% The method uses the polar coordinates form of the integral.
% Gauss-Legendre quadrature with n+1 nodes is used on [0,1], 
% for the integration in the radial variable r.  The number 
%    n = ceil((d1+d2)/2)
%    d1+d2 = degree(Psi(x,y;d1,k1)Psi(x,y;d2,k2))
% The right-hand rectangle rule with 2n+1 nodes is used on
% on [0,2pi] for the angular variable theta.  Because of the
% periodicity in theta of the integrand, this rule is simply
% the trapezoidal rule with 2n+2 nodes and 2n+1 subdivisions
% of the interval [0,2pi]. The choice of n ensures that the
% numerical integration is exact (except for rounding error).
%

n = ceil((d1+d2)/2);
h = 2*pi/(2*n+1);
t = h*[1:2*n+1];
cos_t = cos(t);
sin_t = sin(t);
[r_nodes,r_wgts] = gl_nodes_wts(0,1,n+1);
interm_val = zeros(n+1,1);
for i=1:n+1
    interm_val(i) = sum(gram_fcn(r_nodes(i)*cos_t,r_nodes(i)*sin_t,d1,k1,d2,k2));
end
value = h*sum(r_wgts.*r_nodes.*interm_val);

end % Gram2

% ====================================================================

function ans = gram_fcn(x,y,d1,k1,d2,k2)

% This defines the integrand for the numerical 
% integration over the unit disk as used in Gram2.

ans = delta_ridge_poly(x,y,d1,k1).*delta_ridge_poly(x,y,d2,k2);

end % gram_fcn

% ====================================================================

function ans=ridge_poly(x,y,n,k)

% This evaluates the ridge orthonormal polynomial of
% degree n and index k.

arg = k*pi/(n+1);
c_x = cos(arg); c_y = sin(arg);
z = c_x*x + c_y*y;
ans = (1/sqrt(pi))*cheb2_recurs(z,n);

end % ridge_poly

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

end % cheb2_recurs

% ====================================================================

function ans=delta_ridge_poly(x,y,n,k)

% This evaluates the Laplacian operator applied to
% (1-x^2-y^2)phi(x,y,n,k), with phi(x,y,n,k) the ridge
% orthonormal polynomial of degree n and index k.

arg = k*pi/(n+1);
c_x = cos(arg); c_y = sin(arg);
z = c_x*x + c_y*y;
[U,DU,D2U]=cheb2D_recurs(z,n);
con = 1/sqrt(pi);
ans = -con*(-4*U - 4*z.*DU + (1-(x.^2+y.^2)).*D2U);

end % delta_ridge_poly

% ====================================================================

function [U,DU,D2U]=cheb2D_recurs(x,n)

% This evaluates at x the Chebyshev polynomial of the second
% kind of degree n, along with its first and second derivatives.
% It does so using the triple recursion relation for those
% polynomials.

if n==0
    U = ones(size(x));
    DU = zeros(size(x));
    D2U = zeros(size(x));
    return
end

if n==1
    U = 2*x;
    DU = 2*ones(size(x));
    D2U = zeros(size(x));
    return
end

% n > 1. Use recursion to define U, DU, D2U.

[L1,L2]=size(x);
S = zeros(L1,L2,n+1);       % Dimension
DS = zeros(L1,L2,n+1);      % Dimension
D2S = zeros(L1,L2,n+1);      % Dimension

S(:,:,1) = ones(size(x));
x2 = 2*x;
S(:,:,2) = x2;
for k=2:n
    S(:,:,k+1) = x2.*S(:,:,k) - S(:,:,k-1);
end

DS(:,:,1) = zeros(size(x));
DS(:,:,2) = 2*ones(size(x));
for k=2:n
    DS(:,:,k+1) = 2*S(:,:,k) + x2.*DS(:,:,k) - DS(:,:,k-1);
end

D2S(:,:,1) = zeros(size(x));
D2S(:,:,2) = zeros(size(x));
for k=2:n
    D2S(:,:,k+1) = 4*DS(:,:,k) + x2.*D2S(:,:,k) - D2S(:,:,k-1);
end

U = S(:,:,n+1);
DU = DS(:,:,n+1);
D2U = D2S(:,:,n+1);

end % cheb2D_recurs

%====================================================================

function [nodes,weights]=gl_nodes_wts(a,b,n)

%     [nodes,weights]=gl_nodes_wts(a,b,n)
%
% produces the Gauss-Legendre nodes and weights 
% of order n on the interval [a,b].
%
% a,b     -  integration limits
% n       -  order of formula (default is 20)
% nodes   -  Gauss-Legendre nodes on [a,b]
% weights -  Gauss-Legendre weights on [a,b]
 
u=(1:n-1)./sqrt((2*(1:n-1)).^2-1);
[vc,bp]=eig(diag(u,-1)+diag(u,1));
[bp,k]=sort(diag(bp));
wf=2*vc(1,k)'.^2; 
nodes=(a+b)/2+((b-a)/2)*bp;
weights=(b-a)/2*wf;

end % gl_nodes_wts

% =========================================================

function ans=f_given(x,y,index_f)

switch index_f
case 1
    ans =  10*(1-(x.^2+y.^2));
case 2
    ans =    3 + (x.^2+y.^2);
case 3
    ans =   4 - (1 - (x.^2+y.^2)).^2;
    case 5
        ans = (1-(x.^2+y.^2)).*exp(x.*cos(y));
        case 6
        ans = 1./(1+(x.^2+y.^2));
% case 4
%     ans = u.^2 + delta_ridge_poly(x,y,5,1) - ((1 - (x.^2+y.^2)).*ridge_poly(x,y,5,1)).^2;
% case 5
%     cs = cos(y); sn = sin(y);
%     A = x.*cs;
%     B = exp(A);
%     C = 1 - (x.^2 + y.^2);
%     D = x.*sn;
%     du2x = (-2 - 4*A + C.*cs.^2).*B;
%     du2y = (-2 - A.*C + 4*y.*D + C.*D.^2).*B;
%     delta_u = du2x + du2y;
%     ans = exp(u) - delta_u -exp(B.*C);
end

end % f_given

% =========================================================

% function ans=f_3_given(x,y,u,index_f)
% 
% switch index_f
% case 1
%     ans = -ones(size(u));
% case 2
%     ans =  ones(size(u));
% case 3
%     ans =  2*u;
% case 4
%     ans =  2*u;
% case 5
%     ans =  exp(u);
% end
% 
% end % f_3_given
