%% Brusselator equations on a spherical cap
% This example solves the Brusselator on a spherical cap domain

%exp: Exponential growth

% adjust as appropriate
%addpath('../cp_matrices');
%addpath('../surfaces');

introtime = tic;

format long

%% Construct a grid in the embedding space

R = 1.00;          % cap diameter
gammai = 0.4915;      % initial shape parameter
rhoi = R / gammai;  % initial radius of the sphere
ceni = [0 0 -sqrt(rhoi^2 - R^2)];

gammaf = gammai - 0.04;      % final shape parameter
rhof = R / gammaf;  % final radius of the sphere
cenf = [0 0 -sqrt(rhof^2 - R^2)];

cpfi = @(x,y,z) cpSphereRing(x, y, z, [0 inf], rhoi, ceni);
paramfi = @(N) paramSphereRing(N, [0 inf], rhoi, ceni);

cpff = @(x,y,z) cpSphereRing(x, y, z, [0 inf], rhof, cenf);
paramff = @(N) paramSphereRing(N, [0 inf], rhof, cenf);

gammais = num2str(gammai);
gammafs = num2str(gammaf);

%% time-stepping
Tf = 40000;
dt = 0.1;
numtimesteps = ceil(Tf/dt);
% adjust for integer number of steps
dt = Tf / numtimesteps;
% Steps between growth increments
%gsteps = 50;

%epsilon computation
epsilon = (gammai - gammaf)/Tf;

gamma = gammai - epsilon*dt;
rho = R / gamma;
cen = [0 0 -sqrt(rho^2 - R^2)];

cpf = @(x,y,z) cpSphereRing(x, y, z, [0 inf], rho, cen);
paramf = @(N) paramSphereRing(N, [0 inf], rho, cen);


dx = 0.04;      % grid size ###NOT LOWER THAN 0.025 ON LAPTOP###
dxs = num2str(dx);

% make vectors of x, y, z positions of the grid
x1d = (-1-5*dx:dx:1+5*dx)';
y1d = x1d;
z1d = (-5*dx:dx:1+5*dx)';
nx = length(x1d);
ny = length(y1d);
nz = length(z1d);

% meshgrid is only needed for finding the closest points, not afterwards
[xx, yy, zz] = meshgrid(x1d, y1d, z1d);

% Using "cpbar" [Macdonald, Brandman, Ruuth 2011]:
[cpxi, cpyi, cpzi, disti, bdyi] = cpbar_3d(xx,yy,zz, cpfi); % Initial
cpxi = cpxi(:); cpyi = cpyi(:); cpzi = cpzi(:);

[cpx, cpy, cpz, dist, bdy] = cpbar_3d(xx,yy,zz, cpf); % Initial
cpx = cpx(:); cpy = cpy(:); cpz = cpz(:);

[cpxf, cpyf, cpzf, distf, bdyf] = cpbar_3d(xx,yy,zz, cpff); % Final
cpxf = cpxf(:); cpyf = cpyf(:); cpzf = cpzf(:);


%% Banding: do calculation in a narrow band around the surface
dim = 3;  % dimension
p = 3;    % interpolation order
% "band" is a vector of the indices of the points in the computation
% band.  The formula for bw is found in [Ruuth & Merriman 2008] and
% the 1.0001 is a safety factor.
bw = 1.0001*sqrt((dim-1)*((p+1)/2)^2 + ((1+(p+1)/2)^2));
bandi = find(abs(disti) <= bw*dx);
bandnew = find(abs(dist) <= bw*dx);

% store closest points in the band;
cpxi = cpxi(bandi); cpyi = cpyi(bandi); cpzi = cpzi(bandi);
bdyi = bdyi(bandi);

cpx = cpx(bandnew); cpy = cpy(bandnew); cpz = cpz(bandnew);
bdy = bdy(bandnew);
x = xx(bandnew); y = yy(bandnew); z = zz(bandnew);


%% Toroidal Coordinates
xf = @(t, g, p) R*cos(p).*sinh(t)./(cosh(t) + sqrt(1 - g^2));
yf = @(t, g, p) R*sin(p).*sinh(t)./(cosh(t) + sqrt(1 - g^2));
zf = @(t, g) R*g./(cosh(t) + sqrt(1 - g^2));
eta = log(sqrt((sqrt(cpx.^2 + cpy.^2) + 1).^2 + cpz.^2)./sqrt((sqrt(cpx.^2 + cpy.^2)-1).^2 + cpz.^2));
eta(eta == Inf) = 100;
[phis, rjunk] = cart2pol(cpx, cpy);


%% discrete operators
disp('building laplacian and interp matrices');
L = laplacian_3d_matrix(x1d,y1d,z1d, 2, bandnew, bandnew);
Enew = interp3_matrix(x1d,y1d,z1d, cpx, cpy, cpz, p, bandnew);
E = interp3_matrix(x1d,y1d,z1d, xf(eta, gamma - dt*epsilon, phis), yf(eta, gamma - dt*epsilon, phis), zf(eta, gamma - dt*epsilon), p, bandi);
I = speye(size(Enew));
Ic = I;
Ic(bdy,:) = -I(bdy,:);

% Dirichlet BCs: mirror for ghost points outside of surface edges.
% Comment this out for Neumann BCs.
Enew(bdy,:) = -Enew(bdy,:);

% plotting grid
[xp,yp,zp] = paramf(256);
% Eplot is a matrix which interpolations data onto the plotting grid
Eplot = interp3_matrix(x1d, y1d, z1d, xp(:), yp(:), zp(:), p, bandnew);


%% parameters and functions for Variational Brusselator
% 
%   $$u_t = f(u,g) + nuu*Lap u$$
%   $$v_t = g(u,g) + nuv*Lap u$$
A = 76.51981331;
a = 0.01;  bB = 1.5;  c = 1.8;  d = 0.375;  ddx = 0.005;  ddy = 0.1;
f = @(u,v) ( (bB - d)*u + a^2*A^2*c/d^2*v + bB*d/a/A*u.*u + 2*a*A*c/d*u.*v + c*u.*u.*v);
g = @(u,v) ( -bB*u - a^2*A^2*c/d^2*v - bB*d/a/A*u.*u - 2*a*A*c/d*u.*v - c*u.*u.*v);

%% initial conditions - perturbation from steady state
%pert = 0.5*1500000*besselj(5, 9/R*sqrt(x.^2 + y.^2)).*cos(5*atan2(y, x) + 0.7854) + 0*1000000*rand(size(x));
% From external file
% Remember to match the gamma values with the initial conditions
pert = csvread(['ICs/eigf51_gam_', gammais(3:end), '_dx_', dxs(3:end), '.dat']);% + 0.0*10^8*csvread('eigf03.dat') + 0.0*10^8*csvread('eigf32.dat');

u0 = 0.021*(7.494760741116487)*pert/(max(pert)*7.494760741116487);  
v0 = 0.021*(-0.76388648350)*pert/(max(pert)*7.494760741116487);

u00 = a*A/d*ones(size(bandnew));
v00 = bB*d/(a*A*c)*ones(size(bandnew));
u = u0;  v = v0;



%% Method-of-lines approach
% See [vonGlehn/Macdonald/Maerz 2013]
%lambda = 6*max(nuu,nuv)/(dx^2);
%Au = nuu*(E*L) - lambda*(I-E);
%Av = nuv*(E*L) - lambda*(I-E);

%gmres tolerance and maximal number of iterations
tol = 1e-10;
maxit = 40;


%Linear Matrices - Initial step
Au = I - ddx*dt.*Enew*L + 6*dt/dx^2.*ddx*(I - Enew);
Av = I - ddy*dt.*Enew*L + 6*dt/dx^2.*ddy*(I - Enew);

%Max/min/norm table
maxu = [gammai max(u)];

%%I/O Setup

savedata = 1; %Turn to 1 to save data in csv files

%Creading the data directory
if (savedata == 1)
  timmv = clock;
  timm = [num2str(timmv(1)), '_', num2str(timmv(2)), '_', num2str(timmv(3)), '_', num2str(timmv(4)), '_', num2str(timmv(5)), '_', num2str(floor(timmv(6)))];
  mkdir(['data/', timm]); 
end;

looptime = tic;

%%First Step
uc = E*u;
vc = E*v; %This establishes the boundary conditions on the explicit parts of the computation

[unew, flagu] = gmres(Au,(Ic*uc + dt*(f(uc,vc)) + 2*dt*epsilon/R/sqrt(1 - gamma^2)*(E*cpz).*(uc + u00)), 10, tol, maxit);%, Lu, Uu, u);
[vnew, flagu] = gmres(Av,(Ic*vc + dt*(g(uc,vc)) + 2*dt*epsilon/R/sqrt(1 - gamma^2)*(E*cpz).*(vc + v00)), 10, tol, maxit);%, Lv, Uv, v);

uold = u;
u = unew;
vold = v;
v = vnew;

t = dt;
band = bandi;

for kt = 2:10 %numtimesteps here
  t = kt*dt; %Changed from before
  
  gamma = gammai - kt*dt*epsilon;
  rho = R / gamma;
  cen = [0 0 -sqrt(rho^2 - R^2)];
  
  % Band trickle down
  bandold = band;
  band = bandnew;
    
  % New Surface
  cpf = @(x,y,z) cpSphereRing(x, y, z, [0 inf], rho, cen);
  paramf = @(N) paramSphereRing(N, [0 inf], rho, cen);

  [cpx, cpy, cpz, dist, bdy] = cpbar_3d(xx,yy,zz, cpf); % Initial
  cpx = cpx(:); cpy = cpy(:); cpz = cpz(:);
  
  bandnew = find(abs(dist) <= bw*dx);
  cpx = cpx(bandnew); cpy = cpy(bandnew); cpz = cpz(bandnew);
  bdy = bdy(bandnew);
  x = xx(bandnew); y = yy(bandnew); z = zz(bandnew);
  
  % New toroidal coordinates
  eta = log(sqrt((sqrt(cpx.^2 + cpy.^2) + 1).^2 + cpz.^2)./sqrt((sqrt(cpx.^2 + cpy.^2)-1).^2 + cpz.^2));
  eta(eta == Inf) = 100;
  [phis, rjunk] = cart2pol(cpx, cpy);
  
  % New Operators
  Enew = interp3_matrix(x1d,y1d,z1d, cpx, cpy, cpz, p, bandnew);
  E = interp3_matrix(x1d,y1d,z1d, xf(eta, gamma - dt*epsilon, phis), yf(eta, gamma - dt*epsilon, phis), zf(eta, gamma - dt*epsilon), p, band);
  Eold = interp3_matrix(x1d,y1d,z1d, xf(eta, gamma - 2*dt*epsilon, phis), yf(eta, gamma - 2*dt*epsilon, phis), zf(eta, gamma - 2*dt*epsilon), p, bandold);
  I = speye(size(Enew));
  L = laplacian_3d_matrix(x1d,y1d,z1d, 2, bandnew, bandnew);
  
  Enew(bdy,:) = -Enew(bdy,:);
  Ic = I;
  Ic(bdy,:) = -I(bdy,:);
  
  %disp([size(E) size(cpx) size(u) size(band)])
  %disp([size(Eold) size(cpx) size(uold) size(bandold)])
  
  Au = 3/2*I - ddx*dt.*Enew*L + 6*dt/dx^2.*ddx*(I - Enew);
  Av = 3/2*I - ddy*dt.*Enew*L + 6*dt/dx^2.*ddy*(I - Enew);
  
  % Backtracking u values
  % Matrix sizes don't match because we keep updating them. Use separate
  % variable for computation (uoldc and uc)
  uoldc = Eold*uold;
  uc = E*u;
  voldc = Eold*vold;
  vc = E*v;
  %disp 'bo'
    
  % explicit Euler timestepping
  [unew, flagu] = gmres(Au,(Ic*(2*uc - 1/2*uoldc) + 2*dt*(f(uc,vc)) - dt*(f(uoldc,voldc)) + 2*2*dt*epsilon/R/sqrt(1 - gamma^2)*(E*cpz).*(uc + u00) - 2*dt*epsilon/R/sqrt(1 - gamma^2)*(Eold*cpz).*(uoldc + u00)), 10, tol, maxit);%, Lu, Uu, u);
  [vnew, flagu] = gmres(Av,(Ic*(2*vc - 1/2*voldc) + 2*dt*(g(uc,vc)) - dt*(g(uoldc,voldc)) + 2*2*dt*epsilon/R/sqrt(1 - gamma^2)*(E*cpz).*(vc + v00) - 2*dt*epsilon/R/sqrt(1 - gamma^2)*(Eold*cpz).*(voldc + v00)), 10, tol, maxit);%, Lv, Uv, v);
  
  % storing solutions of the two previous steps
  uold = u;
  u = unew;
  vold = v;
  v = vnew;
  
  if ( (mod(kt,1)==0) )
    %Updating Tables
  %  if kt > 10
        maxu(size(maxu,1)+1,:) = [gamma max(u)];
    %    minu(size(minu,2)+1) = [min(u)];
    %    sumu(size(normu,2)+1) = [sum(u)/size(u,1)];
    %    ttable(size(ttable,2)+1) = [t];
    %
    %    %Updating Plot
    %    figure(2);
    %    Tplot = scatter(ttable, normu);
    %    xlabel('t'); ylabel('max');
  %  end;
    
    %disp([kt t]);
    
    %sphplot = Eplot*u;
    %sphplot = reshape(sphplot, size(xp));
    %set(0, 'CurrentFigure', 1);
    %set(Hplot, 'CData', sphplot);
    %title( ['u at time ' num2str(t) ', kt= ' num2str(kt)] );
    %drawnow;
  end;
  
  if ( (savedata == 1) && (mod(kt, 1)==0) )
    csvwrite(['data/', timm, '/sol_t_', num2str(t), '.dat'], u);
  end;
end

if (savedata == 1)
  dlmwrite(['data/', timm, '/maxu.dat'], maxu, 'delimiter',',','precision',9);
end;

toc(looptime)