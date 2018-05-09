%% Heat equation on a circle
% This example solves the heat equation on a 2D circle, with initial
% conditions u = cos(theta), and exact solution u(t) =
% exp(-t)*cos(theta)
% We also use this file to try various time dependent domiain methods

%torlog_div: Toroidal coordinates growth following logistic growth (increasing xi)
%               Including numerical divergence computation.

% adjust as appropriate
%addpath('../cp_matrices');
%addpath('../surfaces');


%% Construct a grid in the embedding space

dx = 0.04;                   % grid size
R = 1.0/2/pi;
%Rm = 1; %Radius for preconditioner
cen0 = [0 0];

% make vectors of x, y, positions of the grid
x1d = (-5:dx:5)';
y1d = (-2:dx:10)';

nx = length(x1d);
ny = length(y1d);


%% Time-stepping for the heat equation
Tf = 600.0;
%dt = 1*dx;
dt = 0.05;
numtimesteps = ceil(Tf/dt);
% adjust for integer number of steps

%% Domain expansion parameter
rho = 0.01; %Exponential growth factor
zeta = 20; %Relative maximal domain size
%K = rho/numtimesteps; %scaled to kt
%k = rho; %scaled to t
rf = @(t) R*exp(rho*t)/(1 + 1/zeta*(exp(rho*t) - 1));
drf = @(t) R*rho*exp(rho*t)*zeta*(zeta - 1)/(zeta + exp(rho*t) - 1)^2;
xif = @(t) asin(R/rf(t + 0.1)); %xi evolution in time. Matches exponential growth, 0.1 safety time interval
dxif = @(t) -R*drf(t + 0.1)/(rf(t + 0.1)^2*sqrt(rf(t + 0.1)^2 - R^2)); %derivative of xi over time
xi0 = xif(0);
%dxi = 0.00225;
xi = xif(dt);
cen = [0 R*cos(xi)/sin(xi)];

%Logistic growth function
%g = @(t) exp(k*t)/(1 + 1/xi*(exp(k*t) - 1));

%% Find closest points on the surface
% For each point (x,y), we store the closest point on the circle
% (cpx,cpy)

% meshgrid is only needed for finding the closest points, not afterwards
[xx yy] = meshgrid(x1d, y1d);
% function cpCircle for finding the closest points on a circle
[cpx, cpy, dist] = cpCircle(xx,yy,R/sin(xi),cen);
[cpx0, cpy0, dist0] = cpCircle(xx,yy,R/sin(xi0),cen0);
%[cpxm, cpym, distm] = cpCircle(xx,yy,Rm,cenm);
% make into vectors
cpx = cpx(:); cpy = cpy(:);
cpx0 = cpx0(:); cpy0 = cpy0(:);
%cpxm = cpxm(:); cpym = cpy(:);


%% Banding: do calculation in a narrow band around the circle
dim = 2;    % dimension
p = 3;      % interpolation degree
order = 2;  % Laplacian order
% "band" is a vector of the indices of the points in the computation
% band.  The formula for bw is found in [Ruuth & Merriman 2008] and
% the 1.0001 is a safety factor.
bw = 1.0001*sqrt((dim-1)*((p+1)/2)^2 + ((order/2+(p+1)/2)^2));
bandnew = find(abs(dist) <= bw*dx);
band0 = find(abs(dist0) <= bw*dx);
%bandm = find(abs(distm) <= bw*dx);

% store closest points in the band;
cpx = cpx(bandnew); cpy = cpy(bandnew);
cpx0 = cpx0(band0); cpy0 = cpy0(band0);
x = xx(bandnew); y = yy(bandnew);
x0 = xx(band0); y0 = yy(band0);

band = band0;

%% Toroidal Coordinates
xf = @(t, s) real(R*sinh(t)./(cosh(t) - cos(s)));
yf = @(t, s) real(R*sin(s)./(cosh(t) - cos(s)));
eta = sign(cpx).*acosh(0.5*(cpy.*2*pi*cos(xi) + sin(xi))./(cpy*pi));
eta0 = sign(cpx0).*acosh(0.5*(cpy0.*2*pi*cos(xi0) + sin(xi0))./(cpy0*pi));
eta0(isinf(eta0)) = 15;
eta(isinf(eta)) = 15;
dxdxif = @(t, s) real(-R*sinh(t).*sin(s)./(cosh(t) - cos(s)).^2);
dydxif = @(t, s) real(R*(cosh(t).*cos(s) - 1)./(cosh(t) - cos(s)).^2);

%% Function u in the embedding space
% u is a function defined on the grid (eg heat if solving the heat
% equation)
[th, r] = cart2pol(x0,y0);
%u = 0.16*rand(size(x0,1),1) + 1;
u = 1*ones(size(x0,1),1) + 0.16*rand(size(x0,1),1);
v = 0.9*ones(size(x0,1),1);
iu = u;       % store initial value
iv = v;

%% Construct an interpolation matrix for closest point
% This creates a matrix which interpolates data from the grid x1d y1d,
% onto the points cpx cpy.

disp('Constructing interpolation and laplacian matrices');

Enew = interp2_matrix(x1d, y1d, cpx, cpy, p, bandnew);
E = interp2_matrix(x1d, y1d, real(xf(eta, xi0)), real(yf(eta, xi0)), p, band);
I = speye(size(Enew));

%% Velocity field Computation
% We should create a finite difference approximation for the vector field.
% We should use a midpoint method on each computation
%vx = (xf(eta, xif(dt)) - xf(eta, xif(-dt)))/(2*dt);
%vy = (yf(eta, xif(dt)) - yf(eta, xif(-dt)))/(2*dt);

%For precondidioner

%Em = interp2_matrix(x1d, y1d, cpxm, cpym, p, bandm);
%Im = speye(size(Em));

% e.g., closest point extension:
%u = E*u;

%% Create Laplacian matrix and coordinate derivatives for heat equation

L = laplacian_2d_matrix(x1d,y1d, order, bandnew);
[Dxc0, Dyc0] = firstderiv_cen2_2d_matrices(x1d, y1d, band0);
%Lm = laplacian_2d_matrix(x1d,y1d, order, bandm);

%% Schnakenberg Parameters
a = 0.1;
b = 0.9;
d = 0.01;

%% Create combined operator for implicit time-stepping (first step)

Au = I - d*dt.*Enew*L + 4*dt/dx^2.*(I - Enew);
Av = I - dt.*Enew*L + 4*dt/dx^2.*(I - Enew);
%Am = Im - dt.*Em*Lm + 4*dt/dx^2.*(Im - Em);


%% Construct tables for time-position plot
Tplot = []; %Time plot Array
Tplottime = [];

%% plotting grid on circle, using theta as a parameterization

figure(1); clf;
figure(2); clf;
%figure(3); clf;

t = 0;
kt = 0;
thetas = linspace(0,2*pi,100)';
r = R/sin(xi0)*ones(size(thetas));
% plotting grid in Cartesian coords
[xp,yp] = pol2cart(thetas,r);
yp = yp + R*cos(xi0)/sin(xi0);
xp = xp(:); yp = yp(:);
Eplot = interp2_matrix(x1d, y1d, xp, yp, p, band);
  
% plot over computation band
plot2d_compdomain2(u, x0, y0, dx, dx, 1, 0)
title( ['embedded domain: soln at time ' num2str(t) ...
        ', timestep #' num2str(kt)] );
xlabel('x'); ylabel('y');
%hold on
plot(xp,yp,'k-', 'linewidth', 2);
axis equal;  axis tight;

%% gmres tolerance and maximal number of iterations
tol = 1e-12;
maxit = 40;

%% Create directory for files
savedata = 1; %Turn to 1 to save data in csv files

%Creading the data directory
if (savedata == 1)
  timmv = clock;
  timm = [num2str(timmv(1)), '_', num2str(timmv(2)), '_', num2str(timmv(3)), '_', num2str(timmv(4)), '_', num2str(timmv(5)), '_', num2str(floor(timmv(6)))];
  mkdir(['data/', timm]); 
end;

%% First Step
% We first start by sampling the past positions of the closest points
uc = E*u;
vc = E*v;

% We then compute the divergence of the product of the solution with the
% velocity field coordinates.
Vx0 = real(dxdxif(eta0, xif(0)))*real(dxif(0));
Vy0 = real(dydxif(eta0, xif(0)))*real(dxif(0));

%Preconditioner - Initial step
%disp('building preconditioner');
%[Lu,Uu] = ilu(Am,struct('type','ilutp','droptol',5.0e-10));
% TO DO: use correct y values instead of cpy
[unew, flagu] = gmres(Au,uc - dt*E*(Dxc0*(Vx0.*u) + Dyc0*(Vy0.*u)) + dt.*(a - uc + uc.*uc.*vc), 10, tol, maxit);%, Lu, Uu, u);
[vnew, flagv] = gmres(Av,vc - dt*E*(Dxc0*(Vx0.*v) + Dyc0*(Vy0.*v)) + dt.*(b - uc.*uc.*vc), 10, tol, maxit);
uold = u;
u = unew;
vold = v;
v = vnew;

t = dt;
cennew = cen;
cen = cen0;
etanew = eta;
eta = eta0;

%% Create combined operator for implicit time-stepping (subsequent steps)
%A = 3/2*I - dt.*E*L + 4*dt/dx^2.*(I - E);

%Preconditioner - Subsequent steps
%disp('building preconditioner');
%[Lu,Uu] = ilu(A,struct('type','ilutp','droptol',5.0e-10));

for kt = 2:numtimesteps
  % Time and geometry increments
  t = kt*dt;
  
  xi = xif(t);
  
  cenold = cen;
  cen = cennew;
  cennew = [0 R*cos(xi)/sin(xi)];
  
  % Band trickle down
  bandold = band;
  band = bandnew;
  
  % eta trickle down
  etaold = eta;
  eta = etanew;
    
  % New Surface
  [cpx, cpy, dist] = cpCircle(xx,yy,R/sin(xi),cennew);
  cpx = cpx(:); cpy = cpy(:);
  
  bandnew = find(abs(dist) <= bw*dx);
  cpx = cpx(bandnew); cpy = cpy(bandnew);
  x = xx(bandnew); y = yy(bandnew);
  
  etanew = sign(cpx).*acosh(0.5*(cpy.*2*pi*cos(xi) + sin(xi))./(cpy*pi));
  etanew(isinf(etanew)) = 15;
  
  % New Operators
  Enew = interp2_matrix(x1d, y1d, cpx, cpy, p, bandnew);
  E = interp2_matrix(x1d, y1d, real(xf(etanew, xif(t - dt))), real(yf(etanew, xif(t - dt))), p, band);
  Eold = interp2_matrix(x1d, y1d, real(xf(etanew, xif(t - 2*dt))), real(yf(etanew, xif(t - 2*dt))), p, bandold);
  I = speye(size(Enew));
  L = laplacian_2d_matrix(x1d,y1d, order, bandnew);
  [Dxc, Dyc] = firstderiv_cen2_2d_matrices(x1d, y1d, band);
  [Dxcold, Dycold] = firstderiv_cen2_2d_matrices(x1d, y1d, bandold);
  
  %disp([size(E) size(cpx) size(u) size(band)])
  %disp([size(Eold) size(cpx) size(uold) size(bandold)])
  
  Au = 3/2*I - d*dt.*Enew*L + 4*dt/dx^2.*(I - Enew);
  Av = 3/2*I - dt.*Enew*L + 4*dt/dx^2.*(I - Enew);
  
  % Backtracking u values
  % Matrix sizes don't match because we keep updating them. Use separate
  % variable for computation (uoldc and uc)
  uoldc = Eold*uold;
  uc = E*u;
  voldc = Eold*vold;
  vc = E*v;
  %disp 'bo'
  
  % Updating velocity vectors
  Vx = real(dxdxif(eta, xif(t - dt)))*real(dxif(t - dt));
  Vy = real(dydxif(eta, xif(t - dt)))*real(dxif(t - dt));
  Vxold = real(dxdxif(etaold, xif(t - 2*dt)))*real(dxif(t - 2*dt));
  Vyold = real(dydxif(etaold, xif(t - 2*dt)))*real(dxif(t - 2*dt));
    
  % explicit Euler timestepping
  % TO DO: use correct y values instead of cpy
  [unew, flagu] = gmres(Au,(2*uc - 1/2*uoldc) - 2*dt*E*(Dxc*(Vx.*u) + Dyc*(Vy.*u)) + dt*Eold*(Dxcold*(Vxold.*uold) + Dycold*(Vyold.*uold)) + 2*dt.*(a - uc + uc.*uc.*vc) - dt.*(a - uoldc + uoldc.*uoldc.*voldc), 10, tol, maxit);%, Lu, Uu, u);
  [vnew, flagv] = gmres(Av,(2*vc - 1/2*voldc) - 2*dt*E*(Dxc*(Vx.*v) + Dyc*(Vy.*v)) + dt*Eold*(Dxcold*(Vxold.*vold) + Dycold*(Vyold.*vold)) + 2*dt.*(b - uc.*uc.*vc) - dt.*(b - uoldc.*uoldc.*voldc), 10, tol, maxit);
  
  % storing solutions of the two previous steps
  uold = u;
  u = unew;
  vold = v;
  v = vnew;

  if ( (kt < 10) || (mod(kt, ceil(numtimesteps/100)) == 0) || (kt == numtimesteps) )
    % plotting grid on circle, using theta as a parameterization
    thetas = linspace(0,2*pi,100)';
    r = R/sin(xi)*ones(size(thetas));
    % plotting grid in Cartesian coords
    [xp,yp] = pol2cart(thetas,r);
    yp = yp + R*cos(xi)/sin(xi);
    xp = xp(:); yp = yp(:);
    Eplot = interp2_matrix(x1d, y1d, xp, yp, p, bandnew);
      
    % plot over computation band
    plot2d_compdomain2(unew, x, y, dx, dx, 1, 0)
    title( ['embedded domain: soln at time ' num2str(t) ...
            ', timestep #' num2str(kt)] );
    xlabel('x'); ylabel('y');
    %hold on
    plot(xp,yp,'k-', 'linewidth', 2);
    axis equal;  axis tight;
    if mod(kt, ceil(numtimesteps/20)) == 0
        filename = ['data/', timm, '/Plot1_', num2str(kt, '%.4d')];
        print(filename, '-dpng', '-r200');
    end
    
    % plot value on circle
    set(0, 'CurrentFigure', 2);
    clf;
    circplot = Eplot*unew;
    %exactplot = exp(-1/(2*R^2*k))*(exp(exp(-2*k*t)/(2*R^2*k) - k*t))*cos(thetas); %Exact solution for fixed domain
    plot(thetas, circplot);
    title( ['soln at time ' num2str(t) ', on circle'] );
    xlabel('theta'); ylabel('u');
    hold on;
    % plot analytic result
    %plot(thetas, exactplot, 'r--');
    %legend('explicit Euler', 'exact answer', 'Location', 'SouthEast');
    %error_circ_inf = max(abs( exactplot - circplot ));

    %set(0, 'CurrentFigure', 3);
    %clf;
    %plot(thetas, circplot - exactplot);
    %title( ['error at time ' num2str(t) ', on circle'] );
    %xlabel('theta'); ylabel('error');

    %pause(0);
    drawnow();
    if mod(kt, ceil(numtimesteps/20)) == 0
        filename = ['data/', timm, '/Plot2_', num2str(kt, '%.4d')];
        print(filename, '-dpng', '-r200');
    end
  end
  if mod(kt, ceil(numtimesteps/1400)) == 0
    thetas = linspace(0,2*pi,400)';
    r = R/sin(xi)*ones(size(thetas));
    % plotting grid in Cartesian coords
    [xp,yp] = pol2cart(thetas,r);
    yp = yp + R*cos(xi)/sin(xi);
    xp = xp(:); yp = yp(:);
    Eplot = interp2_matrix(x1d, y1d, xp, yp, p, bandnew);
    circplot = Eplot*unew;
    
    %Tplottemp = [t*ones(size(circplot,1),1) R*g(t)*thetas circplot]
    Tplot = [Tplot circplot];
    Tplottime = [Tplottime t];
  end
end

thetas = linspace(0,2*pi,400)';
figure(3); clf;
set(0, 'CurrentFigure', 3);
clf;
fig3 = pcolor(Tplottime, thetas, Tplot); %Look at contourf
set(fig3, 'EdgeColor', 'none');

filename = ['data/', timm, '/Tplot'];
print(filename, '-dpng', '-r1000');

%disp 'relative error'
%max(abs(exactplot - circplot))/max(abs(circplot))

%disp 'absolute error'
%max(abs(exactplot - circplot))