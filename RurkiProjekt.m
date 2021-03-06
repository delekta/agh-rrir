% k -> number of elements, passed as argument
function RurkiProjekt(k) 
    % limits
    x0 = 0;
    xn = 3;
    
    % gravity const
    G = 6.67408 * 10^(-1);
  
    % f(x) 
    f = @(x) (4 * pi * G) .* ((1 < x) * (x <= 2));
    
    % dirichlet shift
    udash = @(x) (-x./3 + 5);
    
    % setting dividing points
    points = linspace(x0, xn, k + 1);
    
    % get base functions and their derivatives
    [ei, eid] = getBaseFunctions(points);
    
    [w, u] = solveDiff(k, f, udash, ei, eid, points);
    
    drawChart(x0, xn, u, w);
end

function [w, u] = solveDiff(k, f, udash, ei, eid, points)
    % derivative of udash
    udashd = @(x) (-1/3);
    
    % number of dividing points
    n = k + 1;
    
    % functions from which we calculate integrals
    B = @(wd,vd, x) (wd(x) .* vd(x));
    Bdash = @(wd,vd, x) (wd(x) .* vd(x)); 
    L = @(v, x) (f(x) .* v(x));
    
    Bmatrix = zeros(n,n);
    Bdashmatrix = zeros(n, 1);
    Lmatrix = zeros(n, 1);
   
    for j = 2:(n - 1)
        for k = 2:(n - 1)
            e1d = eid{j};
            e2d = eid{k};
            
            % function used to integral
            integralFunc = @(x) B(e1d, e2d, x);
            Bmatrix(j,k) = (-(calculateIntegral(integralFunc, points(j-1), points(j)) + calculateIntegral(integralFunc, points(j), points(j+1))));
        end
    end
    
    for i = 2: (n - 1)
        e = ei{i};
        integralFunc2 = @(x) L(e, x);
        Lmatrix(i) = (calculateIntegral(integralFunc2, points(i-1), points(i)) + calculateIntegral(integralFunc2, points(i), points(i+1)));

        ed = eid{i};
        integralFunc3 = @(x) Bdash(udashd, ed, x);
        Bdashmatrix(i) = (calculateIntegral(integralFunc3, points(i-1), points(i)) + calculateIntegral(integralFunc3, points(i), points(i+1)));
    end
  
      % first base function
      e = ei{1};
      ed = eid{1};
      integralFunc = @(x) L(e, x);
      Lmatrix(1) = calculateIntegral(integralFunc, points(1), points(2));
      
      integralFunc2 = @(x) Bdash(udashd, ed, x);
      Bdashmatrix(1) = calculateIntegral(integralFunc2, points(1), points(2));

      integralFunc3 = @(x) B(ed, ed, x);
      Bmatrix(1,1) = (-calculateIntegral(integralFunc3, points(1), points(2)));

      % last base function
      e = ei{n};
      ed = eid{n};
      integralFunc = @(x) L(e, x);
      Lmatrix(n) = calculateIntegral(integralFunc, points(n-1), points(n));
      
      integralFunc2 = @(x) Bdash(udashd, ed, x);
      Bdashmatrix(1) = calculateIntegral(integralFunc2, points(n-1), points(n));
      
      integralFunc3 = @(x) B(ed, ed, x);
      Bmatrix(n,n) = (-calculateIntegral(integralFunc3, points(n-1), points(n)));
      
      % calculatte weights
      wi = Bmatrix \ (Lmatrix - Bdashmatrix);
        
      % scalar producy of weights and base functions
      w = @(x) dot(wi, cellfun((@(a) a(x)), ei));
        
      % shifted w
      u = @(x) w(x) + udash(x);
end

function [functions, derivatives] = getBaseFunctions(xi)
  % number of division points;
  n = length(xi); 
  
  % cell array to store anonymus functions
  functions = cell(n, 1);
  derivatives = cell(n, 1);
  
  % first
  e = @(x) ((xi(2)-x)./(xi(2)-xi(1))) .* ((xi(1)<=x) .* (x<xi(2)));
  ed = @(x) ((-1)./(xi(2)-xi(1))) .* ((xi(1)<=x) .* (x<xi(2)));
  
  functions{1} = e;
  derivatives{1} = ed;
  
  % last
  e = @(x) ((x-xi(n-1))./(xi(n)-xi(n-1))) .* ((xi(n-1)<x) .* (x<=xi(n)));
  ed = @(x) ((1) ./ (xi(n)-xi(n-1))) .* ((xi(n-1)<x) .* (x<=xi(n)));
  
  functions{n} = e;
  derivatives{n} = ed;
  
  for i = 2:(n-1)
    e = @(x)(((x-xi(i-1)) ./ (xi(i)-xi(i-1))) .* ((xi(i-1)<=x) .* (x<xi(i))) + ((xi(i+1)-x)./(xi(i+1)-xi(i))) .* ((xi(i)<=x) .* (x<xi(i+1)))); 
    ed = @(x)(((1) ./ (xi(i)-xi(i-1))) .* ((xi(i-1)<=x) .* (x<xi(i))) + ((-1) ./ (xi(i+1)-xi(i))) .* ((xi(i)<=x) .* (x<xi(i+1))));
    
    functions{i} = e;
    derivatives{i} = ed;
  end
end

% Integral Calculated With Gauss–Legendre Quadrature
function integral = calculateIntegral(f, a, b)
  x1 = -1 / sqrt(3);
  x2 = 1 / sqrt(3);
  
  w1 = 1;
  w2 = 1;

  c1 = (b - a) / 2;
  c2 = (b + a) / 2;

  u1 = c1 * x1 + c2;
  u2 = c1 * x2 + c2;
  
  integral = c1 * (w1 * f(u1) + w2 * f(u2));
end

function drawChart(x0, xn, u, w)
  % Number of points to plot
  argsNum = 1200;

  args = linspace(x0, xn, argsNum);
  
  uValues = zeros(argsNum, 1);
  wValues = zeros(argsNum, 1);

  for i = 1:argsNum
    uValues(i) = u(args(i));
    wValues(i) = w(args(i));
  end
  
  figure();
  hold on;
  plot(args, uValues, 'DisplayName', 'u(x)');
  plot(args, wValues, 'DisplayName', 'w(x)');
  legend;
  xlim([x0, xn])
  hold off;
end