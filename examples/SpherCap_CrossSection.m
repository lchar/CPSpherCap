%% Test file for banding and domain growth

% Domain Parameters and specific functions for initial and final domains

R = 1.00;          % cap diameter
gamma = 0.52;      % initial shape parameter
rho = R / gamma;  % initial radius of the sphere
cen = [0,0];
cen(2) = -sqrt(rho.^2 - R^2); % Initial second coordinate of the circle's centre

gammaf = 0.47;      % final shape parameter
rhof = R / gammaf;  % final radius of the sphere
cenf = [0,0];
cenf(2) = -sqrt(rhof.^2 - R^2); % Final second coordinate of the circle's centre

%Generating the closest point and parameter functions - initial domain
    
cpf = @(x,y) cpArc(x, y, rho, cen, pi/2 - asin(gamma), pi/2 + asin(gamma));
paramf = @(N) paramArc(N, rho, cen, pi/2 - asin(gamma), pi/2 + asin(gamma));
    
%Generating the closest point and parameter functions - final domain
    
cpff = @(x,y) cpArc(x, y, rhof, cenf, pi/2 - asin(gammaf), pi/2 + asin(gammaf));
paramff = @(N) paramArc(N, rhof, cenf, pi/2 - asin(gammaf), pi/2 + asin(gammaf));

% Mesh dimensions

dx = 0.025;      % grid size

% make vectors of x, y, z positions of the grid
x1d = (-2.0:dx:2.0)';
y1d = x1d;
nx = length(x1d);
ny = length(y1d);

% meshgrid is only needed for finding the closest points, not afterwards
[xx yy] = meshgrid(x1d, y1d);

% This function creates closest points at every step from the mesh. Used
% for determining one single band.
[cpxi, cpyi, disti, bdyi] = cpbar_2d(xx,yy, cpf);
cpxi = cpxi(:); cpyi = cpyi(:);

[cpxf, cpyf, distf, bdyf] = cpbar_2d(xx,yy, cpff);
cpxf = cpxf(:); cpyf = cpyf(:);
    

% first-order accurate zero-neumann BC
%[cpx, cpy, cpz, dist, bdy] = cpf(xx,yy,zz);
% Using "cpbar" [Macdonald, Brandman, Ruuth 2011]:

%% Banding: do calculation in a narrow band around the surface
dim = 2;  % dimension
p = 3;    % interpolation order
% "band" is a vector of the indices of the points in the computation
% band.  The formula for bw is found in [Ruuth & Merriman 2008] and
% the 1.0001 is a safety factor.
bw = 1.0001*sqrt((dim-1)*((p+1)/2)^2 + ((1+(p+1)/2)^2));

%Generate the band points - Single uniform band
% ## This method might lead to an incomplete band with missing points in
% ## the middle because this only considers initial and final domain.
% ## Idea: include a binary method using loops.
band = find(abs(disti) <= bw*dx | abs(distf) <= bw*dx);

% keeping only the closest points inside of the band
% ## Remove when using method 1.
cpx = cpxi(band); cpy = cpyi(band);
bdy = bdyi(band);

%Keep only the band points in the mesh
x = xx(band); y = yy(band);

%% Initial discrete operators
disp('building laplacian and interp matrices');
L = laplacian_2d_matrix(x1d,y1d, 2, band, band);
E = interp2_matrix(x1d,y1d, cpx, cpy, p, band);
I = speye(size(E));

% Dirichlet BCs: mirror for ghost points outside of surface edges.
% Comment this out for Neumann BCs.
E(logical(bdy),:) = -E(logical(bdy),:);

%% plotting grid
[xp,yp] = paramf(256);
% Eplot is a matrix which interpolations data onto the plotting grid
Eplot = interp2_matrix(x1d, y1d, xp(:), yp(:), p, band);

%% parameters and functions for Variational Brusselator
% 
%   $$u_t = f(u,g) + nuu*Lap u$$
%   $$v_t = g(u,g) + nuv*Lap u$$
A = 65.0;
a = 0.01;  bB = 1.5;  c = 1.8;  d = 0.375;  ddx = 0.005;  ddy = 0.1;
f = @(u,v) ( (bB - d)*u + a^2*A^2*c/d^2*v + bB*d/a/A*u.*u + 2*a*A*c/d*u.*v + c*u.*u.*v);
g = @(u,v) ( -bB*u - a^2*A^2*c/d^2*v - bB*d/a/A*u.*u - 2*a*A*c/d*u.*v - c*u.*u.*v);

% Time Parameters
epsilon = 10^(-7);
t = 0; % Intial time
T = 0.04/epsilon; % Final time
dt = 0.2 * (1/max(ddx,ddy)) * dx^2;
numtimesteps = ceil(T/dt);
% adjust for integer number of steps
dt = T / numtimesteps;

%% initial conditions - perturbation from steady state
pert =  0.08*rand(size(x));
u0 = pert.*cos(x);  v0 = 0.5*pert.*cos(x);
u = u0;  v = v0;

figure(1);
sphplot = Eplot*u;
sphplot = reshape(sphplot, size(xp));
Hplot = scatter(xp, yp, 25, sphplot, '.');
% caxis([-0.1 0.1]);
title('initial u')
xlabel('x'); ylabel('y');
axis equal
view(0, 90)
%axis off;
shading interp
%camlight left
colorbar

for kt = 1:numtimesteps
    %Generating the closest point and parameter functions
    
    cpf = @(x,y) cpArc(x, y, rho, cen, pi/2 - asin(gamma), pi/2 + asin(gamma));
    paramf = @(N) paramArc(N, rho, cen, pi/2 - asin(gamma), pi/2 + asin(gamma));
    
    % This function creates closest points at every step from the mesh
    % Method 1: Reconstructing closest points from the mesh points inside
    % the band
    %[cpx, cpy, dist, bdy] = cpbar_2d(xx,yy, cpf);
    %cpx = cpx(:); cpy = cpy(:);
    
    % store closest points in the band;
    %cpx = cpx(band); cpy = cpy(band);
    %bdy = bdy(band);

    % Method 2: Generating closest points from previous step's closest
    % points
    % I added a dummy bdy2 because the boundary obtained here will yield
    % error and normally should not have changed, since boundary points are
    % spatially constant
    [cpx, cpy, dist, bdy2] = cpbar_2d(cpx,cpy, cpf);
    cpx = cpx(:); cpy = cpy(:);
        
    %Update Discrete operator E
    L = laplacian_2d_matrix(x1d,y1d, 2, band, band);
    E = interp2_matrix(x1d,y1d, cpx, cpy, p, band);
    I = speye(size(E));
    E(logical(bdy),:) = -E(logical(bdy),:);
    
    %% Ruuth-Merriman
    % Growth version:
    rhsu = ddx*(L*u) + f(u,v) + 2*epsilon/gamma/sqrt(1 - gamma^2)*(cpy).*u;
    rhsv = ddy*(L*v) + g(u,v) + 2*epsilon/gamma/sqrt(1 - gamma^2)*(cpy).*v;
    
    % Without growth term:
    % rhsu = ddx*(L*u) + f(u,v);
    % rhsv = ddy*(L*v) + g(u,v);
    
    unew = u + dt*rhsu;
    vnew = v + dt*rhsv;
    u = E*unew;
    v = E*vnew;
    
    if ( (mod(kt,25)==0) || (kt<=10) || (kt==numtimesteps) )
        %% plotting grid
        [xp,yp] = paramf(256);
        % Eplot is a matrix which interpolations data onto the plotting grid
        Eplot = interp2_matrix(x1d, y1d, xp(:), yp(:), p, band);
        Hplot = scatter(xp, yp, 25, sphplot, '.');
        % disp([kt t]);
        sphplot = Eplot*u;
        sphplot = reshape(sphplot, size(xp));
        set(0, 'CurrentFigure', 1);
        set(Hplot, 'CData', sphplot);
        title( ['u at time ' num2str(t) ', kt= ' num2str(kt)] );
        axis equal
        view(0, 90)
        colorbar
        drawnow;
    end
    
    
    %updating the domain parameters and functions
    gamma = gamma - epsilon*dt;
    rho = R / gamma;
    cen(2) = -sqrt(rho(1).^2 - R^2);
    t = kt*dt;
end