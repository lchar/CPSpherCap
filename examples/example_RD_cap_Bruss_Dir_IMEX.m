%% Reaction-diffusion on a cap
% The cap here is a subset of a sphere: the radius of the cap is
% fixed but as the parameter gamma decreases, the curvature of the
% cap decreases.

R = 1.00;          % cap diameter
gamma = 0.45;      % shape parameter
rho = R / gamma;  % radius of the sphere
cen = [0 0 -sqrt(rho^2 - R^2)];

cpf = @(x,y,z) cpSphereRing(x, y, z, [0 inf], rho, cen);
paramf = @(N) paramSphereRing(N, [0 inf], rho, cen);


loaddata = 1;

if (loaddata == 1)
  dx = 0.05;      % grid size ###NOT LOWER THAN 0.025 ON LAPTOP###

  % make vectors of x, y, z positions of the grid
  x1d = (-2.0:dx:2.0)';
  y1d = x1d;
  z1d = x1d;
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
  
  gc.dim = 3;
  gc.dx = dx;
  gc.x1d = x1d;  gc.y1d = y1d;  gc.z1d = z1d;
  gc.cpfun = cpf;
  gc.band = band;
  gc.x = xx(band);  gc.y = yy(band);  gc.z = zz(band);
  gc.cpx = cpx;  gc.cpy = cpy;  gc.cpz = cpz;
  gc.dist = dist;
  gc.bdy = bdy;
  
  gg = gc;
  %gg = refine_cpgrid_bw(gc,bw);
  
  %x = g.x(g.band); y = g.y(g.band); z = g.z(g.band);


  %% discrete operators
  disp('building laplacian and interp matrices');
  L = laplacian_3d_matrix(gg.x1d,gg.y1d,gg.z1d, 2, gg.band, gg.band);
  E = interp3_matrix(gg.x1d,gg.y1d,gg.z1d, gg.cpx, gg.cpy, gg.cpz, p, gg.band);
  I = speye(size(E));

  % Dirichlet BCs: mirror for ghost points outside of surface edges.
  % Comment this out for Neumann BCs.
  E(gg.bdy,:) = -E(gg.bdy,:);

  %% plotting grid
  [xp,yp,zp] = paramf(256);
  % Eplot is a matrix which interpolations data onto the plotting grid
  Eplot = interp3_matrix(gg.x1d, gg.y1d, gg.z1d, xp(:), yp(:), zp(:), p, gg.band);
end

figure(1); clf;


%% parameters and functions for Variational Brusselator
% 
%   $$u_t = f(u,g) + nuu*Lap u$$
%   $$v_t = g(u,g) + nuv*Lap u$$
A = 76.51981331;
a = 0.01;  bB = 1.5;  c = 1.8;  d = 0.375;  ddx = 0.005;  ddy = 0.1;
f = @(u,v) ( (bB - d)*u + a^2*A^2*c/d^2*v + bB*d/a/A*u.*u + 2*a*A*c/d*u.*v + c*u.*u.*v);
g = @(u,v) ( -bB*u - a^2*A^2*c/d^2*v - bB*d/a/A*u.*u - 2*a*A*c/d*u.*v - c*u.*u.*v);


%% initial conditions - perturbation from steady state
pert = 0.05*exp(-(10*(gg.z-.1)).^2).*cos(5*atan2(gg.y, gg.x)) + 0.08*rand(size(gg.x));
u0 = pert;  v0 = 0.5*pert;
u = u0;  v = v0;

%% time-stepping
Tf = 5000;
dt = 0.1;
%Old time step for explicit method
%dt = 0.1 * (1/max(ddx,ddy)) * gg.dx^2;
numtimesteps = ceil(Tf/dt);
% adjust for integer number of steps
dt = Tf / numtimesteps;


figure(1);
sphplot = Eplot*u;
sphplot = reshape(sphplot, size(xp));
Hplot = surf(xp, yp, zp, sphplot);
title('initial u')
xlabel('x'); ylabel('y'); zlabel('z');
axis equal
view(-10, 60)
%axis off;
shading interp
%camlight left
colorbar

%% Method-of-lines approach
% See [vonGlehn/Macdonald/Maerz 2013]
%lambda = 6*max(nuu,nuv)/(dx^2);
%Au = nuu*(E*L) - lambda*(I-E);
%Av = nuv*(E*L) - lambda*(I-E);

tol = 1e-10;
maxit = 40;


%Linear Matrices
Au = I - ddx*dt.*E*L + 6*dt/dx^2.*(I - E);
Av = I - ddy*dt.*E*L + 6*dt/dx^2.*(I - E);

%Preconditioner
[Lu,Uu] = ilu(I - ddx*dt.*E*L + 6*dt/dx^2.*(I - E),struct('type','ilutp','droptol',1e-4));
[Lv,Uv] = ilu(I - ddy*dt.*E*L + 6*dt/dx^2.*(I - E),struct('type','ilutp','droptol',1e-4));

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
    %HERE! Try gmres with u as initial guess, check restart
    %[unew, flagu] = gmres(I - ddx*dt.*E*L + 6*dt/dx^2.*(I - E),(u + dt*f(u,v)), [], tol, maxit, [], [], u);
    %[vnew, flagv] = gmres(I - ddy*dt.*E*L + 6*dt/dx^2.*(I - E),(v + dt*g(u,v)), [], tol, maxit, [], [], v);
    
    %With flags (for debugging)
    unew = gmres(Au,(u + dt*f(u,v)), 10, tol, maxit, Lu, Uu);
    vnew = gmres(Av,(v + dt*g(u,v)), 10, tol, maxit, Lv, Uv);
    
    %Slash method
    %unew = (I - ddx*dt.*E*L + 6*dt/dx^2.*(I - E))\(u + dt*f(u,v));
    %vnew = (I - ddy*dt.*E*L + 6*dt/dx^2.*(I - E))\(v + dt*g(u,v));
    
    u = unew;
    v = vnew;
    
    t = kt*dt;
    
    disp([kt t]);
    sphplot = Eplot*u;
    sphplot = reshape(sphplot, size(xp));
    set(0, 'CurrentFigure', 1);
    set(Hplot, 'CData', sphplot);
    title( ['u at time ' num2str(t) ', kt= ' num2str(kt)] );
    drawnow;
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
    [unew, flagu] = gmres(Au,(u + dt*f(u,v)), 10, tol, maxit, Lu, Uu);
    [vnew, flagv] = gmres(Av,(v + dt*g(u,v)), 10, tol, maxit, Lv, Uv);
    
    %With flags (for debugging)
    %unew = gmres(I - ddx*dt.*E*L + 6*dt/dx^2.*(I - E),(u + dt*f(u,v)), [], tol, maxit, Lu, Uu, u);
    %vnew = gmres(I - ddy*dt.*E*L + 6*dt/dx^2.*(I - E),(v + dt*g(u,v)), [], tol, maxit, Lv, Uv, v);
    
    %Slash method
    %unew = (I - ddx*dt.*E*L + 6*dt/dx^2.*(I - E))\(u + dt*f(u,v));
    %vnew = (I - ddy*dt.*E*L + 6*dt/dx^2.*(I - E))\(v + dt*g(u,v));
    
    u = unew;
    v = vnew;
    
    t = kt*dt;
  end;
end;

