%% Reaction-diffusion on a cap
% The cap here is a subset of a sphere: the radius of the cap is
% fixed but as the parameter gamma decreases, the curvature of the
% cap decreases.

introtime = tic;

R = 1.00;          % cap diameter
gamma = 0.4915;      % initial shape parameter
rho = R / gamma;  % initial radius of the sphere
cen = [0 0 -sqrt(rho^2 - R^2)];

gammaf = gamma - 0.04;      % final shape parameter
rhof = R / gammaf;  % final radius of the sphere
cenf = [0 0 -sqrt(rhof^2 - R^2)];

cpf = @(x,y,z) cpSphereRing(x, y, z, [0 inf], rho, cen);
paramf = @(N) paramSphereRing(N, [0 inf], rho, cen);

cpff = @(x,y,z) cpSphereRing(x, y, z, [0 inf], rhof, cenf);
paramff = @(N) paramSphereRing(N, [0 inf], rhof, cenf);

gammais = num2str(gamma);
gammafs = num2str(gammaf);

% Separate set of parameters for testing preconditioner range
gammam = gamma - 0.005;      % shape parameter
rhom = R / gammam;  % radius of the sphere
cenm = [0 0 -sqrt(rhom^2 - R^2)];

cpfm = @(x,y,z) cpSphereRing(x, y, z, [0 inf], rhom, cenm);
paramfm = @(N) paramSphereRing(N, [0 inf], rhom, cenm); %Needed for plotting (so not necessary here)

dgamPC = 2*(gamma - gammam); %How much gamma needs to change to switch preconditioners

loaddata = 1;

if (loaddata == 1)
  dx = 0.03;      % grid size ###NOT LOWER THAN 0.025 ON LAPTOP###
  dxs = num2str(dx);

  % make vectors of x, y, z positions of the grid
  x1d = (-1-5*dx:dx:1+5*dx)';
  y1d = x1d;
  z1d = (-5*dx:dx:1+5*dx)';
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
  
  % Same for preconditioner setting
  [cpxm, cpym, cpzm, distm, bdym] = cpbar_3d(xx,yy,zz, cpfm);
  cpxm = cpxm(:); cpym = cpym(:); cpzm = cpzm(:);

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
  
  cpxm = cpxm(band); cpym = cpym(band); cpzm = cpzm(band);
  bdym = bdym(band);
  
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
  
  %Initial Setting for preconditioner
  disp('building laplacian and interp matrices: PC');
  Lm = laplacian_3d_matrix(x1d, y1d, z1d, 2, band, band);
  Em = interp3_matrix(x1d, y1d, z1d, cpxm, cpym, cpzm, p, band);
  Im = speye(size(Em));

  % Dirichlet BCs: mirror for ghost points outside of surface edges.
  % Comment this out for Neumann BCs.
  Em(bdym,:) = -Em(bdym,:);

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
%pert = 0.5*1500000*besselj(5, 9/R*sqrt(x.^2 + y.^2)).*cos(5*atan2(y, x) + 0.7854) + 0*1000000*rand(size(x));
% From external file
% Remember to match the gamma values with the initial conditions
pert = csvread(['ICs/eigf51_gam_0', gammais(3:end), '_0', gammafs(3:end), '_dx_', dxs(3:end), '.dat']);% + 0.0*10^8*csvread('eigf03.dat') + 0.0*10^8*csvread('eigf32.dat');

u0 = 0.021*(7.494760741116487)*pert/(max(pert)*7.494760741116487);  
v0 = 0.021*(-0.76388648350)*pert/(max(pert)*7.494760741116487);

u00 = a*A/d*ones(size(band));
v00 = bB*d/(a*A*c)*ones(size(band));
u = u0;  v = v0;

%Zero initial condition
%u = zeros(size(band));
%v = zeros(size(band));

%% time-stepping
Tf = 40000;
dt = 0.1;
%Old time step formula for explicit method
%dt = 0.1 * (1/max(ddx,ddy)) * gg.dx^2;
numtimesteps = ceil(Tf/dt);
% adjust for integer number of steps
dt = Tf / numtimesteps;
% Steps between growth increments
gsteps = 50;

%epsilon computation
epsilon = (gamma - gammaf)/Tf;

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


%Linear Matrices - Initial step
Au = I - ddx*dt.*E*L + 6*dt/dx^2.*ddx*(I - E);
Av = I - ddy*dt.*E*L + 6*dt/dx^2.*ddy*(I - E);

%Preconditioner - Initial step
disp('building preconditioner');
[Lu,Uu] = ilu(Im - ddx*dt.*Em*Lm + ddx*6*dt/dx^2.*(Im - Em),struct('type','ilutp','droptol',5.0e-4));
[Lv,Uv] = ilu(Im - ddy*dt.*Em*Lm + ddy*6*dt/dx^2.*(Im - Em),struct('type','ilutp','droptol',5.0e-4));
toc (introtime)

%Max/min/norm table
maxu = [gamma max(u)];
%minu = [min(u)];
%sumu = [sum(u)/size(u,1)];
%normu = []; %Use sum instead?
%ttable = [];

%Plot Tables
%figure(2);
%Tplot = scatter(ttable, maxu);
%title('Max Norm')
%xlabel('t'); ylabel('max');

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
[unew, flagu] = gmres(Au,(u + dt*f(u,v) + 2*dt*epsilon/R/sqrt(1 - gamma^2)*cpz.*(u + u00)), 10, tol, maxit, Lu, Uu, u);
[vnew, flagu] = gmres(Av,(v + dt*g(u,v) + 2*dt*epsilon/R/sqrt(1 - gamma^2)*cpz.*(v + v00)), 10, tol, maxit, Lv, Uv, v);

uold = u;
vold = v;

u = unew;
v = vnew;

t = dt;

%Linear Matrices - SBDF2 steps
Au = 3/2*I - ddx*dt.*E*L + 6*dt/dx^2.*ddx*(I - E);
Av = 3/2*I - ddy*dt.*E*L + 6*dt/dx^2.*ddy*(I - E);

%Preconditioner - SBDF2 steps
disp('building preconditioner');
[Lu,Uu] = ilu(3/2*Im - ddx*dt.*Em*Lm + ddx*6*dt/dx^2.*(Im - Em),struct('type','ilutp','droptol',5.0e-4));
[Lv,Uv] = ilu(3/2*Im - ddy*dt.*Em*Lm + ddy*6*dt/dx^2.*(Im - Em),struct('type','ilutp','droptol',5.0e-4));



for kt = 2:numtimesteps/500  %Debugging
  %Growth step
  if mod(kt, gsteps)==0
    %toc(looptime)
    
    %updating the domain parameters and functions
    gamma = gamma - epsilon*gsteps*dt;
    rho = R / gamma;
    cen(3) = -sqrt(rho(1).^2 - R^2);
    
    %Generating the closest point and parameter functions  
    cpf = @(x,y,z) cpSphereRing(x, y, z, [0 inf], rho, cen);
    paramf = @(N) paramSphereRing(N, [0 inf], rho, cen);
    
    % Method 2: Generating closest points from previous step's closest
    % points
    [cpx, cpy, cpz, dist, bdy2] = cpbar_3d(x,y,z, cpf);
    cpx = cpx(:); cpy = cpy(:); cpz = cpz(:);
        
    %Update Discrete operator E
    %disp('Updating laplacian and interp matrices');
    L = laplacian_3d_matrix(x1d,y1d,z1d, 2, band, band);
    E = interp3_matrix(x1d,y1d,z1d, cpx, cpy, cpz, p, band);
    I = speye(size(E));
    E(bdy2,:) = -E(bdy2,:);
    
    %Update Linear Matrices
    Au = 3/2*I - ddx*dt.*E*L + 6*dt/dx^2.*ddx*(I - E);
    Av = 3/2*I - ddy*dt.*E*L + 6*dt/dx^2.*ddy*(I - E);

  end
  if abs(gamma - gammam) > dgamPC/2
    % Update set of parameters for testing preconditioner range
    gammam = gammam - dgamPC;      % shape parameter
    rhom = R / gammam;  % radius of the sphere
    cenm = [0 0 -sqrt(rhom^2 - R^2)];
    
    cpfm = @(x,y,z) cpSphereRing(x, y, z, [0 inf], rhom, cenm);
    paramfm = @(N) paramSphereRing(N, [0 inf], rhom, cenm); %Needed for plotting (so not necessary here)
    
    [cpxm, cpym, cpzm, distm, bdym] = cpbar_3d(xx,yy,zz, cpfm);
    cpxm = cpxm(:); cpym = cpym(:); cpzm = cpzm(:);
    
    cpxm = cpxm(band); cpym = cpym(band); cpzm = cpzm(band);
    bdym = bdym(band);
    
    %Updating Settings for the preconditioner
    disp('Updating laplacian and interp matrices: PC');
    Lm = laplacian_3d_matrix(x1d, y1d, z1d, 2, band, band);
    Em = interp3_matrix(x1d, y1d, z1d, cpxm, cpym, cpzm, p, band);
    Im = speye(size(Em));
    Em(bdym,:) = -Em(bdym,:);
        
    %Update Preconditioner
    disp('building preconditioner');
    [Lu,Uu] = ilu(3/2*Im - ddx*dt.*Em*Lm + ddx*6*dt/dx^2.*(Im - Em),struct('type','ilutp','droptol',5.0e-4));
    [Lv,Uv] = ilu(3/2*Im - ddy*dt.*Em*Lm + ddy*6*dt/dx^2.*(Im - Em),struct('type','ilutp','droptol',5.0e-4));
    %toc
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
  %matrix solver (lsqr, bicg, cgs, gmres)
  %With flags (for debugging)
  %remove flagu and table for convergence info
  %uc = Ic*u;
  %vc = Ic*v;
  %uoldc = Ic*uold;
  %voldc = Ic*vold;%This establishes the boundary conditions on the explicit parts of the computation
  
  [unew, flagu] = gmres(Au,(2*u - 1/2*uold + 2*dt*f(u,v) - dt*f(uold,vold) + 2*2*dt*epsilon/R/sqrt(1 - gamma^2)*cpz.*(u + u00) - 2*dt*epsilon/R/sqrt(1 - gamma^2)*cpz.*(uold + u00)), 10, tol, maxit, Lu, Uu, u);
  [vnew, flagu] = gmres(Av,(2*v - 1/2*vold + 2*dt*g(u,v) - dt*g(uold,vold) + 2*2*dt*epsilon/R/sqrt(1 - gamma^2)*cpz.*(v + v00) - 2*dt*epsilon/R/sqrt(1 - gamma^2)*cpz.*(vold + v00)), 10, tol, maxit, Lv, Uv, v);
  
  %Updating norm table
  %if kt >= 25
  %    normu(size(normu,2)+1) = norm(unew - u);
  %end;
  
  uold = u;
  vold = v;
  
  u = unew;
  v = vnew;
  
  %Advancing time
  t = kt*dt;
  
  %Plotting the results
  if ( (mod(kt,250)==0) || (kt<=10) || (kt==numtimesteps) )
    %Updating Tables
    if kt > 10
        maxu(size(maxu,1)+1,:) = [gamma max(u)];
    %    minu(size(minu,2)+1) = [min(u)];
    %    sumu(size(normu,2)+1) = [sum(u)/size(u,1)];
    %    ttable(size(ttable,2)+1) = [t];
    %
    %    %Updating Plot
    %    figure(2);
    %    Tplot = scatter(ttable, normu);
    %    xlabel('t'); ylabel('max');
    end;
    
    %disp([kt t]);
    
    %sphplot = Eplot*u;
    %sphplot = reshape(sphplot, size(xp));
    %set(0, 'CurrentFigure', 1);
    %set(Hplot, 'CData', sphplot);
    %title( ['u at time ' num2str(t) ', kt= ' num2str(kt)] );
    %drawnow;
  end;
  
  if ( (savedata == 1) && (mod(t, Tf/20)==0) )
    csvwrite(['data/', timm, '/sol_t_', num2str(t), '.dat'], u);
  end;
end;

if (savedata == 1)
  csvwrite(['data/', timm, '/maxu.dat'], maxu);
end;

toc(looptime)