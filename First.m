function [imin emin] = First(x)

i = 0;
imin = 0;
n = 30;
error = 0;
y = 0;
h = 1;
emin = 1;

for i = 1:n
  h = 0.25 * h;
  y = [sin(x -2*h) - 4*sin(x-h) + 6*sin(x) - 4*sin(x+h) + sin(x+2*h)] / h^4;
  error = abs(sin(x) - y);
  
  if (error < emin)
    emin = error;
    imin = i;
  endif
  printf("i: %d h: %e y: %e error: %e\n", i, h, y, error)
endfor
endfunction
