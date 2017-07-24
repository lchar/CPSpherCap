%% Eigenvalue problem on Spherical Cap
% cp_matrices is a folder of useful functions to make implementing the
% closest point method easier. These include closest point extension
% matrices, and differentiation matrices.

% This example computes eigenvalues and eigenfunctions of the
% Laplace--Beltrami operator on a spherical cap.  It also demonstates how
% to impose 2nd-order accurate Neumann and Dirichlet BCs [Macdonald,
% Brandman, Ruuth 2011].


dx = 0.1/2^1;  % grid spacing
dt = 0.05;     % Time step
R=1; %Radius of the circle
gamma = 0.50;      % shape parameter
rho = R / gamma;  % radius of the sphere
cen = [0 0 -sqrt(rho^2 - R^2)];  %Sphere Center

% make vectors of x, y, z positions of the grid
  x1d = (-2.0:dx:2.0)';
  y1d = x1d;
  z1d = x1d;
  nx = length(x1d);
  ny = length(y1d);
  nz = length(z1d);

% meshgrid is only needed for finding the closest points, not afterwards
  [xx yy zz] = meshgrid(x1d, y1d, z1d);

% using the standard CP function, we get a homogeneous Neuamann BC
% [Ruuth & Merriman 2008]
%[cpx,cpy,cpz, dist, bdy] = cpHemisphere(x3d,y3d,z3d R);

% Using "cpbar" [Macdonald, Brandman, Ruuth 2011]:
cpf = @(x,y,z) cpSphereRing(x, y, z, [0 inf], rho, cen);
[cpx,cpy,cpz, dist, bdy] = cpbar_3d(xx,yy,zz, cpf);
paramf = @(N) paramSphereRing(N, [0 inf], rho, cen);

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


%% discrete operators
disp('building laplacian and interp matrices');
L = laplacian_3d_matrix(x1d,y1d,z1d, 2, band,band);
E = interp3_matrix(x1d,y1d,z1d, cpx, cpy, cpz, p, band);

% Dirichlet BCs: mirror for ghost points outside of surface edges.
% Comment this out for Neumann BCs.
E(bdy,:) = -E(bdy,:);

% iCPM matrix
M = lapsharp(L,E, 6*dt/dx^2);


%% plotting grid
[xp,yp,zp] = paramf(256);
%[xp,yp,zp] = paramSphere(32, R);
xp1 = xp(:);  yp1 = yp(:);  zp1 = zp(:);
disp('building plotting matrices');
Eplot = interp3_matrix(x1d,y1d,z1d, xp1,yp1,zp1, p, band);
Eplot0 = interp3_matrix(x1d,y1d,z1d, xp1,yp1,zp1, 0, band);




%% Do some calculation
% Here its eigenvalues/eigenfunctions but could be timestepping
% or whatever using the cp_matrices.

tic
%[V,D] = eigs(-M, 20, 'sm', opts);
[V,D] = eigs(-M, 30, 0.5);
evtime = toc
D = diag(D);
[Lambda,I] = sort(abs(D));
Lambda


tic
figure(2); clf;
for i=1:4
  lambda = Lambda(i+16);
  eigenvec = V(:,I(i+16));

  uplot = Eplot*E*real(eigenvec);
  uplot = reshape(uplot, size(xp));
  subplot(2,2,i);
  surf(xp,yp,zp, uplot);

  xlabel('x'); ylabel('y'); zlabel('z');
  title(['eigenvalue = ' num2str(lambda)]);
  axis equal
  shading interp
  camlight left
end
plottime = toc