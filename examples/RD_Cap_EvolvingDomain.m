%% Reaction-diffusion on a cap
% The cap here is a subset of a sphere: the radius of the cap is
% fixed but as the parameter gamma decreases, the curvature of the
% cap decreases.

R = 1.00;          % cap diameter
gamma = 0.52;      % initial shape parameter
rho = R / gamma;  % initial radius of the sphere
cen = [0 0 -sqrt(rho^2 - R^2)];

epsilon = 10^(-7);

gammaf = 0.47;      % final shape parameter
rhof = R / gammaf;  % final radius of the sphere
cenf = [0 0 -sqrt(rhof^2 - R^2)];

cpf = @(x,y,z) cpSphereRing(x, y, z, [0 inf], rho, cen);
paramf = @(N) paramSphereRing(N, [0 inf], rho, cen);

cpff = @(x,y,z) cpSphereRing(x, y, z, [0 inf], rhof, cenf);
paramff = @(N) paramSphereRing(N, [0 inf], rhof, cenf);


loaddata = 1;

if (loaddata == 1)
  dx = 0.025;      % grid size ###NOT LOWER THAN 0.025 ON LAPTOP###

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
  [cpxi, cpyi, cpzi, disti, bdyi] = cpbar_3d(xx,yy,zz, cpf); % Initial
  cpxi = cpxi(:); cpyi = cpyi(:); cpzi = cpzi(:);
  
  [cpxf, cpyf, cpzf, distf, bdyf] = cpbar_3d(xx,yy,zz, cpff); % Final
  cpxf = cpxf(:); cpyf = cpyf(:); cpzf = cpzf(:);

  %% Banding: do calculation in a narrow band around the surface
  dim = 3;  % dimension
  p = 3;    % interpolation order
  % "band" is a vector of the indices of the points in the computation
  % band.  The formula for bw is found in [Ruuth & Merriman 2008] and
  % the 1.0001 is a safety factor.
  bw = 1.0001*sqrt((dim-1)*((p+1)/2)^2 + ((1+(p+1)/2)^2));
  band = find(abs(disti) <= bw*dx | abs(distf) <= bw*dx);

  % store closest points in the band;
  cpx = cpxi(band); cpy = cpyi(band); cpz = cpzi(band);
  bdy = bdyi(band);
  x = xx(band); y = yy(band); z = zz(band);
  
  %% cprefine step
  % This is a refine step not adapted to the evolving domain yet
  
  %gc.dim = 3;
  %gc.dx = dx;
  %gc.x1d = x1d;  gc.y1d = y1d;  gc.z1d = z1d;
  %gc.cpfun = cpf;
  %gc.band = band;
  %gc.x = xx;  gc.y = yy;  gc.z = zz;
  %gc.cpx = cpx;  gc.cpy = cpz;  gc.cpz = cpz;
  %gc.dist = dist;
  %gc.bdy = bdy;
  
  %gg = refine_cpgrid_bw(gc,bw);
  
  %x = g.x(g.band); y = g.y(g.band); z = g.z(g.band);


  %% discrete operators
  disp('building laplacian and interp matrices');
  L = laplacian_3d_matrix(x1d,y1d,z1d, 2, band, band);
  E = interp3_matrix(x1d,y1d,z1d, cpx, cpy, cpz, p, band);
  I = speye(size(E));

  % Dirichlet BCs: mirror for ghost points outside of surface edges.
  % Comment this out for Neumann BCs.
  E(bdy,:) = -E(bdy,:);

  %% plotting grid
  % Maybe we should use cylindrical coordinates?
  [xp,yp,zp] = paramf(256);
  % Eplot is a matrix which interpolations data onto the plotting grid
  Eplot = interp3_matrix(x1d, y1d, z1d, xp(:), yp(:), zp(:), p, band);
end

figure(1); clf;


%% parameters and functions for Variational Brusselator
% 
%   $$u_t = f(u,g) + nuu*Lap u$$
%   $$v_t = g(u,g) + nuv*Lap u$$
A = 73.0;
a = 0.01;  bB = 1.5;  c = 1.8;  d = 0.375;  ddx = 0.005;  ddy = 0.1;
f = @(u,v) ( (bB - d)*u + a^2*A^2*c/d^2*v + bB*d/a/A*u.*u + 2*a*A*c/d*u.*v + c*u.*u.*v);
g = @(u,v) ( -bB*u - a^2*A^2*c/d^2*v - bB*d/a/A*u.*u - 2*a*A*c/d*u.*v - c*u.*u.*v);


%% initial conditions - perturbation from steady state
pert = 0.05*exp(-(10*(z-.1)).^2).*cos(5*atan2(y, x)) + 0.08*rand(size(x));
u0 = pert;  v0 = 0.5*pert;
u = u0;  v = v0;

%% time-stepping
Tf = 400000;
dt = 0.1 * (1/max(ddx,ddy)) * dx^2;
numtimesteps = ceil(Tf/dt);
% adjust for integer number of steps
dt = Tf / numtimesteps;
% Steps between growth increments
gsteps = 5000;


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

for kt = 1:numtimesteps
  %Growth step
  if mod(kt, gsteps)==0
    %updating the domain parameters and functions
    gamma = gamma - epsilon*dt;
    rho = R / gamma;
    cen(3) = -sqrt(rho(1).^2 - R^2);
    
    %Generating the closest point and parameter functions  
    cpf = @(x,y,z) cpSphereRing(x, y, z, [0 inf], rho, cen);
    paramf = @(N) paramSphereRing(N, [0 inf], rho, cen);
    
    % Method 2: Generating closest points from previous step's closest
    % points
    % I added a dummy bdy2 because the boundary obtained here will yield
    % error and normally should not have changed, since boundary points are
    % spatially constant
    % ##HERE!## Add an if so that the curvature does not change every single step
    [cpx, cpy, cpz, dist, bdy2] = cpbar_3d(cpx,cpy,cpz, cpf);
    cpx = cpx(:); cpy = cpy(:); cpz = cpz(:);
        
    %Update Discrete operator E
    L = laplacian_3d_matrix(x1d,y1d,z1d, 2, band, band);
    E = interp3_matrix(x1d,y1d,z1d, cpx, cpy, cpz, p, band);
    I = speye(size(E));
    E(bdy,:) = -E(bdy,:);
  end
  
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
  % Growth version:
  rhsu = ddx*(L*u) + f(u,v) + 2*epsilon/gamma/sqrt(1 - gamma^2)*(cpz).*u;
  rhsv = ddy*(L*v) + g(u,v) + 2*epsilon/gamma/sqrt(1 - gamma^2)*(cpz).*v;
  unew = u + dt*rhsu;
  vnew = v + dt*rhsv;
  u = E*unew;
  v = E*vnew;

  %Advancing time
  t = kt*dt;
  
  %Plotting the results
  if ( (mod(kt,250)==0) || (kt<=10) || (kt==numtimesteps) )
    %disp([kt t]);
    sphplot = Eplot*u;
    sphplot = reshape(sphplot, size(xp));
    set(0, 'CurrentFigure', 1);
    set(Hplot, 'CData', sphplot);
    title( ['u at time ' num2str(t) ', kt= ' num2str(kt)] );
    drawnow;
  end;
end;

