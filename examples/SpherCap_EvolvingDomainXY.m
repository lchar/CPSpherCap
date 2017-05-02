%% Test file for banding and domain growth
%Issue with Neumann Condition at r=0, solution goes to zero automatically,
%maybe try an averaging method on the extension.

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
    
cpf = @(x,y) cpArc(x, y, rho, cen, pi/2, pi/2 + asin(gamma));
paramf = @(N) paramArc(N, rho, cen, pi/2, pi/2 + asin(gamma));
    
%Generating the closest point and parameter functions - final domain
    
cpff = @(x,y) cpArc(x, y, rhof, cenf, pi/2, pi/2 + asin(gammaf));
paramff = @(N) paramArc(N, rhof, cenf, pi/2, pi/2 + asin(gammaf));

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
    
%% parameters and functions for Variational Brusselator
% 
%   $$u_t = f(u,g) + nuu*Lap u$$
%   $$v_t = g(u,g) + nuv*Lap u$$
A = 85.0; %A = 65 yields a good pattern
a = 0.01;  bB = 1.5;  c = 1.8;  d = 0.375;  ddx = 0.005;  ddy = 0.1;
f = @(X,Y) ( a*A + (-bB - d)*X + c*X.*X.*Y);
g = @(X,Y) ( bB*X - c*X.*X.*Y);

% Time Parameters
epsilon = 10^(-7);
t = 0; % Intial time
T = 0.04/epsilon; % Final time
dt = 0.2 * (1/max(ddx,ddy)) * dx^2/20;
numtimesteps = ceil(T/dt);
% adjust for integer number of steps
dt = T / numtimesteps;

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
L = laplacian_2d_matrix_cylish(x1d, y1d, cpx, cpy, 2, band, band);
EX = interp2_matrix(x1d,y1d, cpx, cpy, p, band);
EY = interp2_matrix(x1d,y1d, cpx, cpy, p, band);
I = speye(size(EX));

% Dirichlet BCs: mirror for ghost points outside of surface edges.
% Comment this out for Neumann BCs.
% We only have Dirichlet for one side (rho = 1). The other boundary is
% Neumann.
%EX(bdy==2,:) = -EX(bdy==2,:);
%EY(bdy==2,:) = -EY(bdy==2,:);
%Non-homogeneous Dirichlet (Does not work)
%EX(bdy==2,:) = 2*a*A/d*logical(EX(bdy==2,:)) - EX(bdy==2,:);
%EY(bdy==2,:) = 2*bB*d/a/A/c*logical(EX(bdy==2,:)) - EY(bdy==2,:);

%% plotting grid
[xp,yp] = paramf(256);
% Eplot is a matrix which interpolations data onto the plotting grid
Eplot = interp2_matrix(x1d, y1d, xp(:), yp(:), p, band);

%% initial conditions - perturbation from steady state
pert =  0.02*rand(size(x));
X0 = pert.*sin(x) + a*A/d;  Y0 = 0.5*pert.*cos(x) + bB*d/a/A/c;
X = X0;  Y = Y0;

%Non-homogeneous Dirichlet BCs
X(bdyi(band)==2) = 2*a*A/d - X(bdyi(band)==2);
Y(bdyi(band)==2) = 2*bB*d/a/A/c - Y(bdyi(band)==2);

figure(1);
sphplot = Eplot*X;
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
    
    cpf = @(x,y) cpArc(x, y, rho, cen, pi/2, pi/2 + asin(gamma));
    paramf = @(N) paramArc(N, rho, cen, pi/2, pi/2 + asin(gamma));
    
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
    [cpx, cpy, dist, bdy] = cpbar_2d(cpx,cpy, cpf);
    cpx = cpx(:); cpy = cpy(:);
        
    %Update Discrete operator E
    L = laplacian_2d_matrix_cylish(x1d, y1d, cpx, cpy, 2, band, band); %##HERE!## Try to compute L directly, like before
    EX = interp2_matrix(x1d,y1d, cpx, cpy, p, band);
    EY = interp2_matrix(x1d,y1d, cpx, cpy, p, band);
    I = speye(size(EX));
    %Homogeneous Dirichlet
    %EX(bdyi(band)==2,:) = - EX(bdyi(band)==2,:);
    %EY(bdyi(band)==2,:) = - EY(bdyi(band)==2,:);
    %Non-homogeneous Dirichlet (Does not work)
    %EX(bdyi(band)==2,:) = 2*a*A/d*logical(EX(bdyi(band)==2,:)) - EX(bdyi(band)==2,:);
    %EY(bdyi(band)==2,:) = 2*bB*d/a/A/c*logical(EX(bdyi(band)==2,:)) - EY(bdyi(band)==2,:);
    
    %% Ruuth-Merriman
    % Growth version:
    rhsX = ddx*(L*X) + f(X,Y) + 2*epsilon/gamma/sqrt(1 - gamma^2)*(cpy).*X;
    rhsY = ddy*(L*Y) + g(X,Y) + 2*epsilon/gamma/sqrt(1 - gamma^2)*(cpy).*Y;
    
    % Without growth term:
    % rhsu = ddx*(L*u) + f(u,v);
    % rhsv = ddy*(L*v) + g(u,v);
    
    Xnew = X + dt*rhsX;
    Ynew = Y + dt*rhsY;
    X = EX*Xnew;
    Y = EY*Ynew;
    
    %Non-homogeneous Dirichlet BCs
    X(bdyi(band)==2) = 2*a*A/d - X(bdyi(band)==2);
    Y(bdyi(band)==2) = 2*bB*d/a/A/c - Y(bdyi(band)==2);
    
    %Neumann Boundary condition at r = 0
    X(xx(band) > -dx/2) = mean(X(abs(xx(band)+dx) < dx/2));
    Y(xx(band) > -dx/2) = mean(Y(abs(xx(band)+dx) < dx/2));
    
    if ( (mod(kt,25)==0) || (kt<=10) || (kt==numtimesteps) )
        %% plotting grid
        [xp,yp] = paramf(256);
        % Eplot is a matrix which interpolations data onto the plotting grid
        Eplot = interp2_matrix(x1d, y1d, xp(:), yp(:), p, band);
        Hplot = scatter(xp, yp, 25, sphplot, '.');
        % disp([kt t]);
        sphplot = Eplot*X;
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