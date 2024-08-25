function E = ECE(Y,x,t,r,varargin)

  % nonparametric density estimation of conditional expectation. Accepts
  % vector input x as well as scalar
  
  % inputs:
  % Y: N by 1 response scalar variable
  % x: N by L matrix, explanatory variable
  % t: M by L to be evaluated, any point sequence inside the span of v, but usually set t=v
  % r: bandwidth
  % varargin: optional input when L > 1, diag_band: =1 if bandwidth is diagonal, neq 1 other
  
N = size(x,1);
M = size(t,1);
L = size(x,2);
E = zeros(M,1);

   if size(t,2) ~= L
       error('regressor dimension mismatch')
   end

   if size(Y,2) ~= 1
       error('Y has to be scalar')
   end
   
%% vector regressor x
if L > 1
    if nargin > 4
        diag_band = varargin{1};
    end
    
    if diag_band == 1
            s = std(x);                 %  1 by L
            h = diag(s)*N^(-r);         %  L by L diagonal bandwidth
        for i = 1:M
        z = (t(i,:)-x)/h;               %  N by L, matrix division
        K = normpdf(z)/h;               %  N by L, normpdf accepts matrix input 
        K = prod(K,2);                  %  N by 1
        E(i) = mean(Y.*K)/mean(K); 
        end
% ---------------------------------------------
    else
    s_2   = x'*x;
    cov_x = s_2./N;
    h     = sqrt(cov_x)*N^(-r);         %  L by L matrix bandwidth
    Sigma = cov_x*N^(-2*r);
        for i = 1:M
            zz = (t(i,:)-x);
            K = mvnpdf(zz,[],Sigma);    %  N by 1 
            % K = prod(K,2);            %  N by 1 
            E(i) = mean(Y.*K)/mean(K); 
        end
    end
    
%% scalar regressor x
else
    s = std(x);
    h = s*N^(-r);                       %  1 by 1 scalar bandwidth
    
        for i = 1:M
        z = (t(i,:)-x)./h;
        K = normpdf(z)./h;      
        K = prod(K,2);                  %  N by 1
        E(i) = mean(Y.*K)/mean(K); 
        end
    
end



end

