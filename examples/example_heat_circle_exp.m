%% Heat equation on a circle
% This example solves the heat equation on a 2D circle, with initial
% conditions u = cos(theta), and exact solution u(t) =
% exp(-t)*cos(theta)
% We also use this file to try various time dependent domiain methods

%exp: Exponential growth

% adjust as appropriate
%addpath('../cp_matrices');
%addpath('../surfaces');


%% Construct a grid in the embedding space

dx = 0.02;                   % grid size
R = 1.0;
%Rm = 1; %Radius for preconditioner
cen = [0 0];
%cenm = [0 0];


% make vectors of x, y, positions of the grid
x1d = (-2.5:dx:2.5)';
y1d = x1d;

nx = length(x1d);
ny = length(y1d);


%% Time-stepping for the heat equation
Tf = 10;
dt = 0.2*dx;
numtimesteps = ceil(Tf/dt);
% adjust for integer number of steps
dt = Tf / numtimesteps;

%% Domain expansion parameter
K = log(2) / (100*numtimesteps); %scaled to kt
k = log(2)/(100*Tf); %scaled to t


%% Find closest points on the surface
% For each point (x,y), we store the closest point on the circle
% (cpx,cpy)

% meshgrid is only needed for finding the closest points, not afterwards
[xx yy] = meshgrid(x1d, y1d);
% function cpCircle for finding the closest points on a circle
[cpx, cpy, dist] = cpCircle(xx,yy,exp(k*dt)*R,cen);
[cpx0, cpy0, dist0] = cpCircle(xx,yy,R,cen);
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


%% Function u in the embedding space
% u is a function defined on the grid (eg heat if solving the heat
% equation)
[th, r] = cart2pol(x0,y0);
u = cos(th);
initialu = u;       % store initial value


%% Construct an interpolation matrix for closest point
% This creates a matrix which interpolates data from the grid x1d y1d,
% onto the points cpx cpy.

disp('Constructing interpolation and laplacian matrices');

Enew = interp2_matrix(x1d, y1d, cpx, cpy, p, bandnew);
E = interp2_matrix(x1d, y1d, exp(-k*dt).*cpx, exp(-k*dt).*cpy, p, band);
I = speye(size(Enew));

%For precondidioner

%Em = interp2_matrix(x1d, y1d, cpxm, cpym, p, bandm);
%Im = speye(size(Em));

% e.g., closest point extension:
%u = E*u;

%% Create Laplacian matrix for heat equation

L = laplacian_2d_matrix(x1d,y1d, order, bandnew);
%Lm = laplacian_2d_matrix(x1d,y1d, order, bandm);

%% Create combined operator for implicit time-stepping (first step)

A = I - dt.*Enew*L + 4*dt/dx^2.*(I - Enew);
%Am = Im - dt.*Em*Lm + 4*dt/dx^2.*(Im - Em);


%% Construct an interpolation matrix for plotting on circle

% plotting grid on circle, using theta as a parameterization
thetas = linspace(0,2*pi,100)';
r = R*ones(size(thetas));
% plotting grid in Cartesian coords
[xp,yp] = pol2cart(thetas,r);
xp = xp(:); yp = yp(:);
Eplot = interp2_matrix(x1d, y1d, xp, yp, p, bandnew);

figure(1); clf;
figure(2); clf;
figure(3); clf;

%% gmres tolerance and maximal number of iterations
tol = 1e-12;
maxit = 40;


%% First Step

uc = E*u;

%Preconditioner - Initial step
%disp('building preconditioner');
%[Lu,Uu] = ilu(Am,struct('type','ilutp','droptol',5.0e-10));

[unew, flagu] = gmres(A,uc - k*dt.*uc, 10, tol, maxit);%, Lu, Uu, u);
uold = u;
u = unew;

t = dt;

%% Create combined operator for implicit time-stepping (subsequent steps)
%A = 3/2*I - dt.*E*L + 4*dt/dx^2.*(I - E);

%Preconditioner - Subsequent steps
%disp('building preconditioner');
%[Lu,Uu] = ilu(A,struct('type','ilutp','droptol',5.0e-10));

for kt = 2:numtimesteps
  
  t = kt*dt; %Changed from before
  
  % Band trickle down
  bandold = band;
  band = bandnew;
    
  % New Surface
  [cpx, cpy, dist] = cpCircle(xx,yy,exp(k*t)*R,cen);
  cpx = cpx(:); cpy = cpy(:);
  
  bandnew = find(abs(dist) <= bw*dx);
  cpx = cpx(bandnew); cpy = cpy(bandnew);
  x = xx(bandnew); y = yy(bandnew);
  
  % New Operators
  Enew = interp2_matrix(x1d, y1d, cpx, cpy, p, bandnew);
  E = interp2_matrix(x1d, y1d, exp(-k*dt).*cpx, exp(-k*dt).*cpy, p, band);
  Eold = interp2_matrix(x1d, y1d, exp(-2*k*dt).*cpx, exp(-2*k*dt).*cpy, p, bandold);
  I = speye(size(Enew));
  L = laplacian_2d_matrix(x1d,y1d, order, bandnew);
  
  %disp([size(E) size(cpx) size(u) size(band)])
  %disp([size(Eold) size(cpx) size(uold) size(bandold)])
  
  A = 3/2*I - dt.*Enew*L + 4*dt/dx^2.*(I - Enew);
  
  % Backtracking u values
  % Matrix sizes don't match because we keep updating them. Use separate
  % variable for computation (uoldc and uc)
  uoldc = Eold*uold;
  uc = E*u;
  %disp 'bo'
    
  % explicit Euler timestepping
  [unew, flagu] = gmres(A,(2*uc - 1/2*uoldc) - 2*k*dt.*uc + k*dt.*uoldc, 10, tol, maxit);%, Lu, Uu, u);

  % storing solutions of the two previous steps
  uold = u;
  u = unew;

  if ( (kt < 10) || (mod(kt,1/10/dt) == 0) || (kt == numtimesteps) )
    % plotting grid on circle, using theta as a parameterization
    thetas = linspace(0,2*pi,100)';
    r = R*exp(k*t)*ones(size(thetas));
    % plotting grid in Cartesian coords
    [xp,yp] = pol2cart(thetas,r);
    xp = xp(:); yp = yp(:);
    Eplot = interp2_matrix(x1d, y1d, xp, yp, p, bandnew);
      
    % plot over computation band
    plot2d_compdomain2(unew, x, y, dx, dx, 1, 0)
    title( ['embedded domain: soln at time ' num2str(t) ...
            ', timestep #' num2str(kt)] );
    xlabel('x'); ylabel('y');
    %hold on
    plot(xp,yp,'k-', 'linewidth', 2);
    %axis equal;  axis tight;

    % plot value on circle
    set(0, 'CurrentFigure', 2);
    clf;
    circplot = Eplot*unew;
    exactplot = exp(-k*t)*(exp((exp(-2*k*t) - 1)/(2*R^2*k)))*cos(thetas); %Exact solution for fixed domain
    plot(thetas, circplot);
    title( ['soln at time ' num2str(t) ', on circle'] );
    xlabel('theta'); ylabel('u');
    hold on;
    % plot analytic result
    plot(thetas, exactplot, 'r--');
    legend('explicit Euler', 'exact answer', 'Location', 'SouthEast');
    error_circ_inf = max(abs( exactplot - circplot ));

    set(0, 'CurrentFigure', 3);
    clf;
    plot(thetas, circplot - exactplot);
    title( ['error at time ' num2str(t) ', on circle'] );
    xlabel('theta'); ylabel('error');

    %pause(0);
    drawnow();
  end
end

format long
disp 'relative error'
max(abs(exactplot - circplot))/max(abs(circplot))

disp 'absolute error'
max(abs(exactplot - circplot))