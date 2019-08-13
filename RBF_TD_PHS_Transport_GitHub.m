function W = RBF_TD_PHS_Transport_GitHub(x,t,xc,tc,alpha,tscale,m,d)

% Input parameters
%   x,t     Column vectors with the x- and t-coordinates within the stencil
%   xc,tc   Stencil center location
%   alpha   Local wave speed
%   m       Order of PHS
%   d       Polynomial order to attach
%   tscale  Scaling of the time dimension
% Output parameters
%     W     RBF-FD weights for u_t+alpha*u_x

%Shift nodes
x = x - xc; 
t = t - tc;
n = length(x);

fi  =@(x,t) (x.^2+(tscale*t).^2).^(m/2);     %PHS; Radial function
Lfi =@(x,t) -m*(tscale^2*t+alpha*x).*(x.^2+(tscale*t).^2).^((m-2)/2); %Operator acting on PHS

A11 = fi(x-x',t-t');                % RBF matrix within stencil

L1 = Lfi(x,t); 

if d == -1                     % Special case; no polynomial terms,
    A = A11;  L = L1;    
else                        % Include polynomials and matching constraints
  
    X    =  x(:,ones(1,d+1));  X(:,1) = 1;  X = cumprod( X,2);
    Y    =  t(:,ones(1,d+1));  Y(:,1) = 1;  Y = cumprod( Y,2);
    np   = (d+1)*(d+2)/2;            % Total number of polynomial terms
    A12  = zeros(n,np);  col = 1;   % Assemble polynomial matrix block
    L2 = zeros(np,1); % Create the remaining part of the RHS vector L  
    for k = 0:d                 % Create A12 block
        A12(:,col:col+k) = X(:,k+1:-1:1).*Y(:,1:k+1);
        col = col+k+1;
    end  
    
    if(d>=1) L2(2)=alpha; L2(3)=1; end
    
    L = [L1;L2];     % Assemble the linear system to be solved
    
    A = [A11,A12;A12',zeros(np)];
    
end

w = A\L;
W = w(1:n);