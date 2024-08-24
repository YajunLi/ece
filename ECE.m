function E = ECE(Y,v,t,r)

  % Y: N by 1 response scalar variable
  % v: N by 1 matrix, explanatory vector
  % t: M by 1 to be evaluated
  % r: param to bandwidth
  % t could be any point sequence inside the span of v, but usually set t=v
  
N = size(v,1);
M = size(t,1);
s = std(v)';   % column vector 
h = s*N^(-r);  % vector bandwidth
E = zeros(M,1);

 i = 1;
    while i <= M
        z = (t(i,:)-v)./h;
        K = normpdf(z)./h;
        K = prodc(K');
        E(i) = mean(Y.*K)/mean(K); 
        i = i + 1;
    end

end


