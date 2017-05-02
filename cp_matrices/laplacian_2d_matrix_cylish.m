function L = laplacian_2d_matrix(x,y, cpx, cpy, order, band1, band2, use_ndgrid, use_loop)
%LAPLACIAN_2D_MATRIX_CYLISH  Build a 2D discrete proto-cylindrical Laplacian
%
%   DEPRECATED?: you may want laplacian_matrix()
%
%   L = laplacian_2d_matrix(x, y, order, band)
%      'L' is a discrete laplacian over a grid.
%      'order' can be 2 or 4.
%      x,y are 1D vectors which form a meshgrid, 'band' is a
%      subset of this meshgrid given as linear indices.
%      TODO: this will have a dirichlet boundary condition at the
%      edge of the band, at least for order 2.
%
%   L = laplacian_3d_matrix(x, y, order, band1, band2)
%      dual-banded version of above, a rectangular matrix.
%      TODO: explain this
%
%   To use ndgrid ordering pass "true" as an extra argument.
%
%   Pass "true" as a further argument to use the (older, slower)
%   looping-based code.
%
%   Currently if dx is within 10*macheps of dy, this assumes dx==dy,
%   otherwise it uses the more general formula.  The logic is to have
%   marginally less rounding error in the common dx == dy case: is it
%   worth it?
%
%   Does no error checking on the equispaced nature of x,y.
%
%   We try here to reproduce cylindrical coordinates. Check leeds.ac.uk
%   page on "Finite differences in polar coordinates" on how to handle the
%   radial direction and the singularity at the centre of the cap.

  
  m = 5; % Emulates the effect of the m = 5 mode
  
  if (nargin <= 8)
    use_loop = false;
  end
  if (nargin <= 7)
    use_ndgrid = false;
  end
  if (nargin <= 6)
    band2 = band1;
  end

  % input checking
  [temp1, temp2] = size(x);
  if ~(  (ndims(x) == 2) && (temp1 == 1 || temp2 == 1)  )
    error('x must be a vector, not e.g., meshgrid output');
  end
  [temp1, temp2] = size(y);
  if ~(  (ndims(y) == 2) && (temp1 == 1 || temp2 == 1)  )
    error('y must be a vector, not e.g., meshgrid output');
  end

  dx = x(2)-x(1);
  dy = y(2)-y(1);
  if assertAlmostEqual(dx, dy, 10*eps)
    dxequal = 1;
  else
    dxequal = 0;
  end

  if (order == 2)
    if dxequal
      weights = [-4 1 1 1 1] / dx^2;
      weights2 = [0 0 0 1/2 -1/2] / dx; % 1/r d/dr part of cylindrical laplacian
      weights3 = [-m^2 0 0 0 0]; % Emulation of the angle derivative for m = 5 mode
    else
      weights = [ -2/dx^2 - 2/dy^2   [1 1]/dx^2   [1 1]/dy^2 ];
    end

    PTS = [ 0   0; ...
            1   0; ...
           -1   0; ...
            0   1; ...
            0  -1];
  elseif (order == 4)
    if dxequal
      weights = [-5.0 ...
                 (-1.0/12.0)  (4.0/3.0)  (4.0/3.0)  (-1.0/12.0) ...
                 (-1.0/12.0)  (4.0/3.0)  (4.0/3.0)  (-1.0/12.0) ...
                ] / dx^2;
    else
      weights = [ -5/(2*dx^2) - 5/(2*dy^2) ...
                [-1/12  4/3  4/3  -1/12] / dx^2  ...
                [-1/12  4/3  4/3  -1/12] / dy^2 ];
    end
    PTS = [ 0   0; ...
           -2   0; ...
           -1   0; ...
            1   0; ...
            2   0; ...
            0  -2; ...
            0  -1; ...
            0   1; ...
            0   2];
  else
    error(['order ' num2str(order) ' not implemented']);
  end

  % Bypassing the Helper file
  
  Nx = length(x);
  Ny = length(y);

  StencilSize = length(weights);
  
  Li = repmat((1:length(band1))', 1, StencilSize);
  Lj = zeros(size(Li));
  Ls = zeros(size(Li));
  
  [i,j] = ind2sub([Nx,Ny], band1);
  
  for c = 1:StencilSize
    ii = i + PTS(c,1);
    jj = j + PTS(c,2);
    Ls(:,c) = weights(c) + weights2(c)./abs(cpx(floor(band1/Nx+1))) + weights3(c)./(cpx(floor(band1/Nx+1))).^2;
    Lj(:,c) = sub2ind([Nx,Ny],ii,jj);
  end
  
  % Singularity at r = 0
  % ###HERE!### Issues with r=0 point. The two sides of the solutions seem
  % disconnected, maybe reduce problem to only one half
  Ls(floor(band1/Nx+1) == ceil(Nx/2), :) = repmat([-6 2 2 1 1] / dx^2, length(Ls(floor(band1/Nx+1) == ceil(Nx/2), :)), 1);
    
  L = sparse(Li(:), Lj(:), Ls(:), length(band1), Nx*Ny);

  % TODO: these sorts of checks could move to the ops and bands replacement
  % If we're using careful banding a la iCPM2009 then as a sanity
  % check all of the columns outside of band2 should be zero.
  % if (~isempty(ICPM2009BANDINGCHECKS)) && (ICPM2009BANDINGCHECKS)
  %   Lout = L(:, setdiff(1:(Nx*Ny),band2));
  %   if (nnz(Lout) > 0)
  %     nnz(Lout)
  %     error('Lost some non-zero coefficients (from outside the outerband)');
  %   end
  % end

  % remove columns not in band2 (the outerband)
  L = L(:, band2);

  %Ltime = toc;


