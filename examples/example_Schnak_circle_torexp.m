%% Heat equation on a circle
% This example solves the heat equation on a 2D circle, with initial
% conditions u = cos(theta), and exact solution u(t) =
% exp(-t)*cos(theta)
% We also use this file to try various time dependent domiain methods

%torexp: Toroidal coordinates growth following exponential growth (increasing xi)

% adjust as appropriate
%addpath('../cp_matrices');
%addpath('../surfaces');


%% Construct a grid in the embedding space

dx = 0.02;                   % grid size
R = 1.0/2/pi;
%Rm = 1; %Radius for preconditioner
cen0 = [0 0];

% make vectors of x, y, positions of the grid
x1d = (-5:dx:5)';
y1d = (-2:dx:10)';

nx = length(x1d);
ny = length(y1d);


%% Time-stepping for the heat equation
Tf = 300;
%dt = 1*dx;
dt = 0.1;
numtimesteps = ceil(Tf/dt);
% adjust for integer number of steps

%% Domain expansion parameter
rho = 0.01; %Exponential growth factor
%K = rho/numtimesteps; %scaled to kt
%k = rho; %scaled to t
xif = @(t) asin(exp(-rho*(t +0.1))); %xi evolution in time. Matches exponential growth
dxif = @(t) -rho*exp(-rho*(t+0.1))/sqrt(1 - (exp(-rho*(t+0.1)))^2); %derivative of xi over time
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
xf = @(t, s) R*sinh(t)./(cosh(t) - cos(s));
yf = @(t, s) R*sin(s)./(cosh(t) - cos(s));
eta = sign(cpx).*acosh(0.5*(cpy.*2*pi*cos(xi) + sin(xi))./(cpy*pi));

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

%For precondidioner

%Em = interp2_matrix(x1d, y1d, cpxm, cpym, p, bandm);
%Im = speye(size(Em));

% e.g., closest point extension:
%u = E*u;

%% Create Laplacian matrix for heat equation

L = laplacian_2d_matrix(x1d,y1d, order, bandnew);
%Lm = laplacian_2d_matrix(x1d,y1d, order, bandm);

%% Schnakenberg Parameters
a = 0.1;
b = 0.9;
d = 0.01;

%% Create combined operator for implicit time-stepping (first step)

Au = I - d*dt.*Enew*L + 4*dt/dx^2.*(I - Enew);
Av = I - dt.*Enew*L + 4*dt/dx^2.*(I - Enew);
%Am = Im - dt.*Em*Lm + 4*dt/dx^2.*(Im - Em);


%% Construct an interpolation matrix for plotting on circle

% plotting grid on circle, using theta as a parameterization
thetas = linspace(0,2*pi,100)';
r = R*ones(size(thetas));
% plotting grid in Cartesian coords
[xp,yp] = pol2cart(thetas,r);
xp = xp(:); yp = yp(:);
Eplot = interp2_matrix(x1d, y1d, xp, yp, p, bandnew);

Tplot = []; %Time plot Array
Tplottime = [];

figure(1); clf;
figure(2); clf;
%figure(3); clf;

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

uc = E*u;
vc = E*v;

%Preconditioner - Initial step
%disp('building preconditioner');
%[Lu,Uu] = ilu(Am,struct('type','ilutp','droptol',5.0e-10));
% TO DO: use correct y values instead of cpy
[unew, flagu] = gmres(Au,uc + dxif(dt)*dt/R*cpy.*uc + dt.*(a - uc + uc.*uc.*vc), 10, tol, maxit);%, Lu, Uu, u);
[vnew, flagv] = gmres(Av,vc + dxif(dt)*dt/R*cpy.*vc + dt.*(b - uc.*uc.*vc), 10, tol, maxit);
uold = u;
u = unew;
vold = v;
v = vnew;

t = dt;
cennew = cen;
cen = cen0;

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
    
  % New Surface
  [cpx, cpy, dist] = cpCircle(xx,yy,R/sin(xi),cennew);
  cpx = cpx(:); cpy = cpy(:);
  
  bandnew = find(abs(dist) <= bw*dx);
  cpx = cpx(bandnew); cpy = cpy(bandnew);
  x = xx(bandnew); y = yy(bandnew);
  
  eta = sign(cpx).*acosh(0.5*(cpy.*2*pi*cos(xi) + sin(xi))./(cpy*pi));
  
  % New Operators
  Enew = interp2_matrix(x1d, y1d, cpx, cpy, p, bandnew);
  E = interp2_matrix(x1d, y1d, real(xf(eta, xif(t - dt))), real(yf(eta, xif(t - dt))), p, band);
  Eold = interp2_matrix(x1d, y1d, real(xf(eta, xif(t - 2*dt))), real(yf(eta, xif(t - 2*dt))), p, bandold);
  I = speye(size(Enew));
  L = laplacian_2d_matrix(x1d,y1d, order, bandnew);
  
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
    
  % explicit Euler timestepping
  % TO DO: use correct y values instead of cpy
  [unew, flagu] = gmres(Au,(2*uc - 1/2*uoldc) + 2*dxif(t - dt)*dt/R*cpy.*uc - dxif(t - 2*dt)*dt/R*cpy.*uoldc + 2*dt.*(a - uc + uc.*uc.*vc) - dt.*(a - uoldc + uoldc.*uoldc.*voldc), 10, tol, maxit);%, Lu, Uu, u);
  [vnew, flagv] = gmres(Av,(2*vc - 1/2*voldc) + 2*dxif(t - dt)*dt/R*cpy.*vc - dxif(t - 2*dt)*dt/R*cpy.*voldc + 2*dt.*(b - uc.*uc.*vc) - dt.*(b - uoldc.*uoldc.*voldc), 10, tol, maxit);
  
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

%disp 'relative error'
%max(abs(exactplot - circplot))/max(abs(circplot))

%disp 'absolute error'
%max(abs(exactplot - circplot))