function driver_cg_source(N0,dt,kappa,xmax,xmaxout,test)
  format compact
  diary driver_cg_source.dia
  rem = mod(log2(N0-1),1); % remainder of k=log2(N0-1), i.e., if k is integer
  if (rem ~= 0.0)
    error('N0 must be of form 2^k+1 such as 33, 65, 129, 257, 513, 1025, etc.');
  end;
  if (test ~= 'y' & test ~= 'n')
    error('test must be either ''y'' or ''n''.');
  end;
  S.test = test;
  S.N0 = N0;
  S.dt = dt;
  N = N0^2;
  S.N = N;

  % ... set up problem parameters for domain
  %     (xmin,xmax) x (ymin,ymax) with specific symmetric choices
  %     xmin=-xmax, ymin=-ymax, and xmax=ymax,
  %     and then L = xmax - xmin = length of domain in one dimension:
  S.xmax = xmax;
  S.xmin = -S.xmax;
  S.ymax = S.xmax;
  S.ymin = -S.ymax;
  S.L = S.xmax - S.xmin;
  S.xmaxout = xmaxout;
  % ... other physical parameters:
  S.tau = 8.0;
  S.D = 10.0;
  S.kappa = kappa;
  S.t_fin = 24.0; % final time

  % ... compute the number of time steps needed to solve from t=0 to t=t_fin:
  Nt = ceil(S.t_fin/dt); % ceil rounds up, if needed
  S.Nt = Nt;
  % ... compute uniform mesh spacing and set up spatial mesh:
  h = S.L / (N0-1);
  S.h = h;
  S.x = [S.xmin : h : S.xmax];
  S.y = S.x;
  [X, Y] = ndgrid(S.x,S.y);
  A = setupA(N0);
  A = speye(N) + (S.D*(dt/h^2)) * A;
  tol = 1.0e-9;
  S.tol = tol;
  maxit = 999;
  S.maxit = maxit;

  U = zeros(N0,N0,6); %Array to store desired solutions
  u_n = zeros(N0,N0);          % solution at initial time
  u_n = u_n(:); % reshape for linear solves throughout

  fprintf('driver_cg_source:\n');
  fprintf('kappa =%14.6e D =%14.6e tau =%14.6e\n', S.kappa, S.D, S.tau);
  fprintf('xmax =%14.6e xmaxout =%14.6e t_fin =%14.6e\n', ...
          S.xmax, S.xmaxout, S.t_fin);
  fprintf('N0 =%6d N =%12d h =%14.6e\n', S.N0, S.N, S.h);
  fprintf('dt =%14.6e Nt =%14.6e\n', S.dt, S.Nt);
  fprintf('tol =%14.6e maxit =%5d\n', S.tol, S.maxit);
  
  fprintf(['     n  t_n  it cumit   ', ...
           ' min(U(:))   max(U(:))  mass(U(:))    enorminf\n']);

  tout = 4.0; % first output time
  dtout = 4.0; % output time step thereafter
  vtout = [tout : dtout : S.t_fin]; % vector of output times
  ntout = 1;
  tout = vtout(ntout);
  cum_iter = 0; % cumulative number of iterations
  plot_max = 0.0;
  for n = 0 : Nt
      t = n*dt; % current or old time
      tnew = (n+1)*dt; % = t + dt = new time, at which PDE is evaluated

      F = frhs(X,Y,tnew,S);
      b = u_n + dt * F(:);

      % ... CG with same u_n vector as initial guess and solution:
      [u_n,flag,relres,iter,resvec] = pcg(A,b,tol,maxit,[],[],u_n);
      cum_iter = cum_iter + iter;

      % ... output every dtout times, up to 6 times total:
      if ( (tnew>=tout-1.0e-14) && (ntout<=6) )
          U_n = reshape(u_n, [N0 N0]); % reshape solution only if output needed
          U(:,:,ntout) = U_n(:,:);
          Utrue = utrue(X,Y,tnew,S);
          E = U_n - Utrue;
          enorminf = max(abs(E(:)));
          min_u = min(u_n(:));
          max_u = max(u_n(:));
          V = U_n(1:end-1,1:end-1) ...
            + U_n(2:end  ,1:end-1) ...
            + U_n(1:end-1,2:end  ) ...
            + U_n(2:end  ,2:end  );
          mass_t = (h/2)^2 * sum(V(:));
          fprintf('%6d%5.1f%4d%6d%13.4e%12.4e%12.4e%12.4e\n', ...
                  n, tnew, iter, cum_iter, min_u, max_u, mass_t, enorminf);
          if (ntout > 2) % consider max_u only for times t > 8 hours
            if (max_u > plot_max)
              plot_max = max_u;
            end;
          end;
          ntout = ntout + 1;
          if (ntout<=6), tout = vtout(ntout); end;
      end
  end

  % ... overwrite X, Y, U array by their own center portions, namely such that
  %     its domain in x and y is (-xmaxout,xmaxout)x(-xmaxout,xmaxout):
  %     Example: numerical xmax=800, output xmaxout=50;
  %     then numerics are on (-xmax,xmax)x(-xmax,xmax) domain, and I want to
  %     cut out the small center on (-xmaxout,xmaxout)x(-xmaxout,xmaxout).
  %     The problem is that I need to act on the indices into the
  %     arrays X, Y, and U. So, I need to explicitly add and subtract
  %     from the center index. For clarity, I obtain sizes manually.
  %     Assumptions: size(X)=size(Y)=size(U(:,:,k)) for all k
  %     and that X, Y, U(:,:,k) are square. They are actually of size N0-by-N0.
  %     m is an odd integer in the following, like N0=2^k+1.
disp('[min(x), max(x)], size(X) before selecting center area:')
[min(X(:)), max(X(:))]
size(X)
  [m, notused] = size(X);
  mc = (m+1)/2; % = center index
  r = xmax/xmaxout; % = ratio = 800/50 = 16 for example
  ms = floor(mc-(m-1)/(2*r)); % = start index
  me = ceil(mc+(m-1)/(2*r)); % = end index
  %     Ideally, want me-ms+1=2^k+1 for positive integer k;
  %     so, me-mc=2^k desired
  len = me-mc;
  loglen = ceil(log2(len));
  ms = mc - 2^loglen;
  me = mc + 2^loglen;
  X = X(ms:me,ms:me);
  Y = Y(ms:me,ms:me);
  U = U(ms:me,ms:me,:);
disp('[min(x), max(x)], size(X) after selecting center area:')
[min(X(:)), max(X(:))]
size(X)

  % ... compute output mesh spacing so that only a 33x33 mesh is plotted:
  %     Example: m=65, N0out=33, then stride s=(m-1)/(N0out-1)=64/32=2
  [m, notused] = size(X);
  N0out = 33;
  if (N0out < m)
    s = (m-1) / (N0out-1); % = stride
  else
    s = 1;
  end;
  Xout = X(1:s:end,1:s:end);
  Yout = Y(1:s:end,1:s:end);
  % make plot_max integer, if it is larger than 1, and 1 otherwise:
  % Notice this assumes a certain reasonable scale of our solutions.
  if (plot_max >= 1.0)
    plot_max = round(plot_max);
  else
    plot_max = 1.0;
  end;
  fprintf('Mesh plots with plot_max =%14.6e s =%4d\n', plot_max, s);
  % plot numerical solution:
  figure('DefaultAxesFontSize',10);
  for k = 1 : 6
      t = vtout(k);
      subplot(2,3,k)
      H = mesh(Xout,Yout,U(1:s:end,1:s:end,k));
      min_u = min(min(U(:,:,k)));
      max_u = max(max(U(:,:,k)));
      title_str = sprintf('t = %g', t);
      title(title_str);
      xlabel('x');
      ylabel('y');
      zlabel(sprintf('%.1e<=u_h<=%.1e',min_u,max_u));
      zlim([-inf plot_max]);
  end
  % orient tall;
  % saveas(H,'Solution.png')
  print -dpng solution.png

  % compute and plot numerical error:
  figure('DefaultAxesFontSize',10);
  for k = 1 : 6
      t = vtout(k);
      Utrue = utrue(X,Y,t,S);
      E = U(:,:,k) - Utrue;
      subplot(2,3,k)
      H = mesh(Xout,Yout,E(1:s:end,1:s:end));
      title_str = sprintf('t = %g',t);
      title(title_str);
      xlabel('x');
      ylabel('y');
      zlabel('u-u_h');
  end
  % orient tall;
  % saveas(H,'Error.png')
  print -dpng error.png

  % compute L^inf norm of error and print:
  %enorminf = max(abs(E(:)));
  %fprintf('N = %5d\n', N);
  %fprintf('tol = %10.1e, maxit = %d\n', tol, maxit);
  %fprintf('flag = %1d, iter = %d, relres = %24.16e\n', flag, iter, relres);
  %fprintf('h                  = %24.16e\n', h);
  %fprintf('h^2                = %24.16e\n', h^2);
  %fprintf('enorminf           = %24.16e\n', enorminf);
  %fprintf('C = enorminf / h^2 = %24.16e\n', (enorminf/h^2));
  %fprintf('wall clock time    = %10.2f seconds\n', timesec);

  diary off;

function A = setupA(N0)
  I = speye(N0);
  s = [-1*ones(1,N0-1) 2*ones(1,N0) -1*ones(1,N0-1)]';
  s( 1 ) = -2;
  s(end) = -2;
  i = [1:N0-1  1:N0  2:N0  ]';
  j = [2:N0    1:N0  1:N0-1]';
  T = sparse(i,j,s);
  A = kron(I,T) + kron(T,I);

function F = frhs(X,Y,t,S)
  if S.test == 'y'
    F = ((2*t/S.tau^2)*exp(-(t/S.tau)^2)) ...
          * ((cos((pi/S.L)*X)).^2.*(cos((pi/S.L)*Y)).^2) ...
      - (S.D * (1-exp(-(t/S.tau)^2)) * (-2*(pi/S.L)^2)) ...
          * ( cos((2*pi/S.L)*X).*(cos((pi/S.L)*Y)).^2 ... 
            + (cos((pi/S.L)*X)).^2.*cos((2*pi/S.L)*Y) );
  else 
    F = zeros(size(X));
      if (t <= 8)
       F((end+1)/2,(end+1)/2) = S.kappa / (S.h^2);
      end;
  end;

function U = utrue(X,Y,t,S)
  if S.test == 'y'
    U = (1-exp(-(t/S.tau)^2)) * ((cos((pi/S.L)*X)).^2.*(cos((pi/S.L)*Y)).^2);
  else
    U = zeros(size(X));
  end;

