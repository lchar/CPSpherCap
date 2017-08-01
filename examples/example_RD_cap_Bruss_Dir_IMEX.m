%% Reaction-diffusion on a cap
% The cap here is a subset of a sphere: the radius of the cap is
% fixed but as the parameter gamma decreases, the curvature of the
% cap decreases.

tic

R = 1.00;          % cap diameter
gamma = 0.45;      % shape parameter
rho = R / gamma;  % radius of the sphere
cen = [0 0 -sqrt(rho^2 - R^2)];

cpf = @(x,y,z) cpSphereRing(x, y, z, [0 inf], rho, cen);
paramf = @(N) paramSphereRing(N, [0 inf], rho, cen); %Needed for plotting


loaddata = 1;

if (loaddata == 1)
  dx = 0.05;      % grid size ###NOT LOWER THAN 0.025 ON LAPTOP (perhaps not true anymore)###

  % make vectors of x, y, z positions of the grid
  x1d = (-1.5:dx:1.5)';
  y1d = x1d;
  z1d = (-0.5:dx:1.5)';
  nx = length(x1d);
  ny = length(y1d);
  nz = length(z1d);

  % meshgrid is only needed for finding the closest points, not afterwards
  [xx yy zz] = meshgrid(x1d, y1d, z1d);

  % first-order accurate zero-neumann BC
  %[cpx, cpy, cpz, dist, bdy] = cpf(xx,yy,zz);
  % Using "cpbar" [Macdonald, Brandman, Ruuth 2011]:
  [cpx, cpy, cpz, dist, bdy] = cpbar_3d(xx,yy,zz, cpf);

  cpx = cpx(:); cpy = cpy(:); cpz = cpz(:);

  %% Banding: do calculation in a narrow band around the surface
  dim = 3;  % dimension
  p = 3;    % interpolation order
  % "band" is a vector of the indices of the points in the computation
  % band.  The formula for bw is found in [Ruuth & Merriman 2008] and
  % the 1.0001 is a safety factor.
  bw = 1.0001*sqrt((dim-1)*((p+1)/2)^2 + ((1+(p+1)/2)^2));
  band = find(abs(dist) <= bw*dx);
  %band = ( 1:length(xx(:)) )';

  % store closest points in the band;
  cpx = cpx(band); cpy = cpy(band); cpz = cpz(band);
  x = xx(band); y = yy(band); z = zz(band);
  bdy = bdy(band);
  
  %cprefine step
  
  dim = 3;

  %% discrete operators
  disp('building laplacian and interp matrices');
  L = laplacian_3d_matrix(x1d, y1d, z1d, 2, band, band);
  E = interp3_matrix(x1d, y1d, z1d, cpx, cpy, cpz, p, band);
  I = speye(size(E));

  % Dirichlet BCs: mirror for ghost points outside of surface edges.
  % Comment this out for Neumann BCs.
  E(bdy,:) = -E(bdy,:);
 
  %% plotting grid
  [xp,yp,zp] = paramf(256);
  % Eplot is a matrix which interpolations data onto the plotting grid
  Eplot = interp3_matrix(x1d, y1d, z1d, xp(:), yp(:), zp(:), p, band);
end

%figure(1); clf;


%% parameters and functions for Variational Brusselator
% 
%   $$u_t = f(u,g) + nuu*Lap u$$
%   $$v_t = g(u,g) + nuv*Lap u$$
A = 76.51981331;
a = 0.01;  bB = 1.5;  c = 1.8;  d = 0.375;  ddx = 0.005;  ddy = 0.1;
f = @(u,v) ( (bB - d)*u + a^2*A^2*c/d^2*v + bB*d/a/A*u.*u + 2*a*A*c/d*u.*v + c*u.*u.*v);
g = @(u,v) ( -bB*u - a^2*A^2*c/d^2*v - bB*d/a/A*u.*u - 2*a*A*c/d*u.*v - c*u.*u.*v);


%% initial conditions - perturbation from steady state
pert = 0.5*1000000*besselj(5, 9/R*sqrt(x.^2 + y.^2)).*cos(5*atan2(y, x) + 0.7854) + 1*1000000*rand(size(x));
% From external file
%pert = 0.3*10^8*csvread('eigf51.dat') + 0.3*10^8*csvread('eigf03.dat') + 0.3*10^8*csvread('eigf32.dat');

u0 = 1.0*(7.494760741116487 * 0.009046337986/(9.046337990295378/0.00001))*pert;  
v0 = 1.0*(-0.76388648350 * 0.009046337986/(9.046337990295378/0.00001))*pert;
u = u0;  v = v0;

%% time-stepping
Tf = 5;
dt = 0.1;
%Old time step for explicit method
%dt = 0.1 * (1/max(ddx,ddy)) * gg.dx^2;
numtimesteps = ceil(Tf/dt);
% adjust for integer number of steps
dt = Tf / numtimesteps;


%figure(1);
%sphplot = Eplot*u;
%sphplot = reshape(sphplot, size(xp));
%Hplot = surf(xp, yp, zp, sphplot);
%title('initial u')
%xlabel('x'); ylabel('y'); zlabel('z');
%axis equal
%view(-10, 60)
%axis off;
%shading interp
%camlight left
%colorbar

%% Method-of-lines approach
% See [vonGlehn/Macdonald/Maerz 2013]
%lambda = 6*max(nuu,nuv)/(dx^2);
%Au = nuu*(E*L) - lambda*(I-E);
%Av = nuv*(E*L) - lambda*(I-E);

%gmres tolerance and maximal number of iterations
tol = 1e-10;
maxit = 40;


%Linear Matrices
Au = I - ddx*dt.*E*L + 6*dt/dx^2.*ddx*(I - E);
Av = I - ddy*dt.*E*L + 6*dt/dx^2.*ddy*(I - E);

%Preconditioner
%tic
[Lu,Uu] = ilu(I - ddx*dt.*E*L + ddx*6*dt/dx^2.*(I - E),struct('type','ilutp','droptol',1.25e-2));
[Lv,Uv] = ilu(I - ddy*dt.*E*L + ddy*6*dt/dx^2.*(I - E),struct('type','ilutp','droptol',1.25e-2));
%toc

%Max/min/norm table
%maxu = [max(u)];
%minu = [min(u)];
%sumu = [sum(u)/size(u,1)];
%normu = []; %Use sum instead?
%ttable = [];

%Plot Tables
%figure(2);
%Tplot = scatter(ttable, maxu);
%title('Max Norm')
%xlabel('t'); ylabel('max');

for kt = 1:numtimesteps
  %% MOL: explicit Euler timestepping
  %unew = u + dt*( E*f(u,v) + Au*u );
  %vnew = v + dt*( E*g(u,v) + Av*v );
  %u = unew;
  %v = vnew;
  %% MOL: without precomputing matrices
  %rhsu = nuu*(L*u) + f(u,v);
  %rhsv = nuv*(L*v) + g(u,v);
  %unew = u + dt*( E*rhsu - lambda*(u - E*u) );
  %vnew = v + dt*( E*rhsv - lambda*(v - E*v) );
  %u = unew;
  %v = vnew;

  %% Ruuth-Merriman
  
  if ( (mod(kt,25)==0) || (kt<=10) || (kt==numtimesteps) )
    %Old explicit method: Will need to use this when evolving the domain
    %with a smaller time step.
    %rhsu = ddx*(L*u) + f(u,v);
    %rhsv = ddy*(L*v) + g(u,v);
    %unew = u + dt*rhsu;
    %vnew = v + dt*rhsv;
    %u = E*unew;
    %v = E*vnew;
    
    %Preconditioner
    %[Lu,Uu] = ilu(I - ddx*dt.*E*L + 6*dt/dx^2.*(I - E),struct('type','ilutp','droptol',1e-10));
    %[Lv,Uv] = ilu(I - ddy*dt.*E*L + 6*dt/dx^2.*(I - E),struct('type','ilutp','droptol',1e-10));
    
    %matrix solver (lsqr, bicg, cgs, gmres)
    %With flags (for debugging)
    %remove flagu and table for convergence info
    [unew, flagu] = gmres(Au,(u + dt*f(u,v)), 10, tol, maxit, Lu, Uu, u);
    [vnew, flagu] = gmres(Av,(v + dt*g(u,v)), 10, tol, maxit, Lv, Uv, v);

    %Updating norm table
    %if kt >= 25
    %    normu(size(normu,2)+1) = norm(unew - u);
    %end;
    
    u = unew;
    v = vnew;
    
    t = kt*dt;
    
    %disp([kt t]);
    %sphplot = Eplot*u;
    %sphplot = reshape(sphplot, size(xp));
    %set(0, 'CurrentFigure', 1);
    %set(Hplot, 'CData', sphplot);
    %title( ['u at time ' num2str(t) ', kt= ' num2str(kt)] );
    %drawnow;
    
    %Updating Tables
    %if kt >= 25
    %    maxu(size(maxu,2)+1) = [max(u)];
    %    minu(size(minu,2)+1) = [min(u)];
    %    sumu(size(normu,2)+1) = [sum(u)/size(u,1)];
    %    ttable(size(ttable,2)+1) = [t];
    %
    %    %Updating Plot
    %    figure(2);
    %    Tplot = scatter(ttable, normu);
    %    xlabel('t'); ylabel('max');
    %end;
  else
    %Old explicit method
    %rhsu = ddx*(L*u) + f(u,v);
    %rhsv = ddy*(L*v) + g(u,v);
    %unew = u + dt*rhsu;
    %vnew = v + dt*rhsv;
    %u = E*unew;
    %v = E*vnew;
    
    %IMEX method: Laplacian and penalty handled implicitly and reaction
    %handled explicietly
    %rhsu = ddx*(L*unew) - 6/dx^2.*(I - E)*unew + f(u,v);
    %rhsv = ddy*(L*vnew) - 6/dx^2.*(I - E)*vnew + g(u,v);
    %unew = u + dt*rhsu;
    %vnew = v + dt*rhsv;
    %u = E*unew;
    %v = E*vnew;
    
    %(I - ddx*dt.*E*L + 6*dt/dx^2.*(I - E)))*unew = u + dt*f(u,v)
    
    %Preconditioner
    %[Lu,Uu] = ilu(I - ddx*dt.*E*L + 6*dt/dx^2.*(I - E),struct('type','ilutp','droptol',1e-6));
    %[Lv,Uv] = ilu(I - ddy*dt.*E*L + 6*dt/dx^2.*(I - E),struct('type','ilutp','droptol',1e-6));

    %inv() doesn't cut it here, we will need to implement an iterative
    %matrix solver (lsqr, bicg, cgs, gmres)
    [unew, flagu] = gmres(Au,(u + dt*f(u,v)), 10, tol, maxit, Lu, Uu, u);
    [vnew, flagv] = gmres(Av,(v + dt*g(u,v)), 10, tol, maxit, Lv, Uv, v);
    
    u = unew;
    v = vnew;
    
    t = kt*dt;

  end;
end;

toc