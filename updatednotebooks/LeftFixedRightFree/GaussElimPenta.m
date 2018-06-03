function [X] = GaussElimPenta (E, A, D, C, F, B)
n = length(B);
X = zeros(n,1);
for k = 1:n-2
  xmult = A(k)/D(k);
  D(k+1) = D(k+1) - xmult*C(k);
  C(k+1) = C(k+1) - xmult*F(k);
  B(k+1) = B(k+1) - xmult*B(k);
  
  xmult = E(k)/D(k);
  A(k+1) = A(k+1) - xmult*C(k);
  D(k+2) = D(k+2) - xmult*F(k);
  B(k+2) = B(k+2) - xmult*B(k);
endfor

xmult = A(n-1)/D(n-1);
D(n) = D(n) - xmult*C(n-1);
B(n) = B(n) - xmult*B(n-1);

X(n) = B(n)/D(n);
X(n-1) = (B(n-1)-C(n-1)*X(n))/D(n-1);

for i = n-2:-1:1
  X(i) = (B(i) - C(i)*X(i+1)-F(i)*X(i+2))/D(i);
endfor
endfunction

