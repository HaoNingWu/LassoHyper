function [err_I,errest_I,err2_I,LV,LCalfa,Ltime,Calfa,f_name,Finterp]=CX_Interpolation(n,f_idx,sigma)
% Based on M-function "CX_Interpolation" of Stefano De Marchi and Marco Vianello, we 
% write this function for filtered hyperinterpolation in a cube

%---- Grid of target points ------------
m= 5;
targets=linspace(-1,1,m);     
[xix,xiy,xiz]=meshgrid(targets);
Ptar=[xix(:) xiy(:) xiz(:)];
%--------------------------------------

t=cputime;
c=cos([0:n]*pi/n);     % Chebsyhev points
[X,Y,Z]=ndgrid(c);     % 3d grid of Chebsyhev points in the cube
Sc=[X(:), Y(:), Z(:)]; % this is the Chebyshev grid arranged as an array with 3 coordinates

% Here we extract the sub-vectors E=EVEN and O=ODD of Chebsyhev-Lobatto points 
E=c(1:2:length(c));
O=c(2:2:length(c));

%----------------------------------------
% generation of the sub-grids 
%
% We have 4 possible choices
% 1) EEE OOO
% 2) EOO OEE
% 3) OEE EOO
% 4) EOE OEO
%----------------------------------------
[X1,Y1,Z1]=meshgrid(O,E,E);
[X2,Y2,Z2]=meshgrid(E,O,O);

EEE=[X1(:),Y1(:),Z1(:)];
OOO=[X2(:),Y2(:),Z2(:)];
U=[EEE;OOO];
%our_npunti(n-n0+1)=size(U,1);

[MU,IU]=setdiff(Sc,U,'rows'); % indeces of the points of product Chebyshev-Lobatto grid
                              % which arae not Xu points on the Cube
                         
% The function on the Xu points and zero elsewhere 
LSc=length(Sc(:,1));
f=zeros(LSc,1);
JU=setdiff(1:LSc,IU);
V=Sc(JU,:);  % Xu points on the cube in 3D
% [f(JU),fname]=fCube(V(:,1),V(:,2),V(:,3));


x=V(:,1); y=V(:,2);z=V(:,3);
switch f_idx
    case 1
        f(JU) = exp(-(x.^2+y.^2+z.^2));
        f_name='$\exp\left(-(x^2+y^2+z^2)\right)$';
    case 2
        f(JU)=1./(1+16*(x.^2+y.^2+z.^2));
        f_name='$1/(1+16(x^2+y^2+z^2))$';
    case 3
        f(JU) = sqrt((x.^2+y.^2+z.^2)).^3;
        f_name='$(x^2+y^2+z^2)^{3/2}$';
    case 4
        f(JU) = exp(abs(x+y+z));
        f_name='$\exp(|x+y+z|)$';
    case 5
        f(JU) = exp(-1./(x.^2+y.^2+z.^2));
        f_name='$\exp(-1/(x^2+y^2+z^2))$';
end
        
        
% f(find(f~=0)) = noisegen(f(find(f~=0)),10);
f(find(f~=0)) = f(find(f~=0))+sigma.*randn(size(find(f~=0)));



%----------------------
% compute the weights
%----------------------

[F,W]= CX_weights(n,f,c,V,JU);

%-------------------------------
% Compute the coefficients Calfa
% ------------------------------
[Calfa,M,P,Calfa1]=CX_Calfa(F,n);  


errest_I=sum(abs(Calfa1));
% ----------------------------------
% Construction of the interpolant
%----------------------------------
ind=M(P,:); % matrix of indices

%---- Interpolation on a grid of target points ------------

scaling=[1,ones(1,n)*sqrt(2)]';  % scaling for Chebyshev polynomials

m=length(Calfa);
D=diag(scaling);

TX=D*cheb(Ptar(:,1)',n);
TY=D*cheb(Ptar(:,2)',n);
TZ=D*cheb(Ptar(:,3)',n);

%Palfa=(D*cheb(Ptar(ind(1:m,1),1)',n)).*(D*cheb(Ptar(ind(1:m,2),2)',n)).*(D*cheb(Ptar(ind(1:m,3),3)',n));

Palfa=TX(ind(:,1),:).*TY(ind(:,2),:).*TZ(ind(:,3),:);

for i = 1:m
   Calfa(i) = Filter(n,ind(i,3))*Calfa(i);
end
    


L=Palfa'*Calfa; % the vector containing the values of the interpolant 

% f_targets=fCube(Ptar(:,1),Ptar(:,2),Ptar(:,3));
x=Ptar(:,1); y=Ptar(:,2);z=Ptar(:,3);


Finterp=scatteredInterpolant(x,y,z,L);


switch f_idx
    case 1
        f_targets = exp(-(x.^2+y.^2+z.^2));
    case 2
        f_targets=1./(1+16*(x.^2+y.^2+z.^2));
    case 3
        f_targets = sqrt((x.^2+y.^2+z.^2)).^3;
    case 4
        f_targets = exp(abs(x+y+z));
    case 5
        f_targets = exp(-1./(x.^2+y.^2+z.^2));
end


f_mean=sum(f_targets)/length(Ptar(:,1));
err_I=norm(L-f_targets, inf);%/max(abs(f_targets-f_mean));  % rel error in infinity norm

[err2_I]= sqrt(CX_integral(4,(L-f_targets).^2));


LCalfa=length(Calfa);
LV=length(V(:,1));
Ltime=cputime-t;




