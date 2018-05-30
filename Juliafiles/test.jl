L = 10; # length = 10 m
b = 0.1; # width = 10 cm
d = 0.05; # height = 5 cm
E = 2*10^11; # Young�s modulus for steel = 200 GPa = 2x10^11 Pa
I = b*d^3/12; # second moment of inertia
rho = 7850; # mass density of steel = 7850 kg/m^3
g = 9.81; # acceleration due to gravity = 9.81 m/s^2
w = rho*b*d*g; # weight of the beam per unit length (will be our f)

block_rho = 750; # the mass density of wood
block_length = 1; # in meters, the length of the block is 1m
block_height = 1; # in meters, the height of the block is .2m, or 20cm

# to conform to the width of the width of the board so the whole block is in contact with the board,
# (using the mass density of steel now)
# the block's width will simply be the boards width, b.  Computing the mass, we get:
block_mass = block_rho*block_height*block_length*b;

# Alternatively, we can compute the mass per unit length
block_mass_per_unit = block_rho*block_height*g*b;

n = 10; # number of subintervals on [0, L]
h = L/n; # discretization spacing
N = n + 1; # number of unknowns to be solved for
A = spzeros(N,N); # generating a sparse matrix

# Define the RHS of the system
f = -h^4/(E*I) * w * ones(N, 1);
f[1] = 0;
f[N] = 0;
# Creating diagonals of the matrix
for i=3:(N - 2)
    A[i,i] = 6;
    A[i,i-1] = A[i,i+1] = -4;
    A[i,i-2] = A[i,i+2] = 1; 
end

# Left end
A[1,1] = 1;
A[2,2] = 7;
A[1,2] = 0;
A[1,3] = 0;
A[2,1] = 0;
A[3,1] = 0;
A[2,3] = -4;
A[2,4] = 1;
# Right end
A[N,N] = 1;
A[N-1,N-1] =  7;
A[N-1,N] = 0;
A[N-2,N] = 0;
A[N, N-1] = 0;
A[N, N-2] = 0;
A[N-1,N-2] = -4;
A[N-1,N-3] = 1;

# Solve for y
y = A\f;
x = ones(N,1);

for i=1:N
    x[i] = i * h;
end
y_exact = -b*d*rho*g/(24*E*I)*x.^2.*(L - x).^2;
ErrMax = maximum(abs.(y-y_exact))
display(ErrMax)