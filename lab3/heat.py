from fenics import *
from mshr import *

T = 10000.0            # final time
num_steps = 100     # number of time steps
dt = T / num_steps # time step size
alpha = 3          # parameter alpha
beta = 1.2         # parameter beta

# Create mesh and define function space
c1 = Cylinder(Point(0, 0,0), Point(0,0,3),8,8,50)
c2 = Cylinder(Point(0, 0,0),Point(0,0,3),5, 5,50)
domain = c1-c2
mesh = generate_mesh(domain,45)
V = FunctionSpace(mesh, 'P', 1)

# Define boundary condition
u_D = Expression('100/(10*(x[0]-6.5*cos(pi*t/1000))*(x[0]-6.5*cos(pi*t/1000))+10*(x[1]-6.5*sin(pi*t/1000))*(x[1]-6.5*sin(pi*t/1000))+10.0*(x[2]-1.5)*(x[2]-1.5)+1)',element=V.ufl_element(),t=0)

def boundary(x, on_boundary):
    return on_boundary

bc = DirichletBC(V, u_D, boundary)

# Define initial value
u_n = interpolate(u_D, V)
#u_n = project(u_D, V)

# Define variational problem
u = TrialFunction(V)
v = TestFunction(V)
f = Constant(beta - 2 - 2*alpha)

F = u*v*dx + dt*dot(grad(u), grad(v))*dx - (u_n + dt*f)*v*dx
a, L = lhs(F), rhs(F)

res_file = File('heat/solution.pvd')

# Time-stepping
u = Function(V)
t = 0
for n in range(num_steps):

    # Update current time
    t += dt
    u_D.t = t

    # Compute solution
    solve(a == L, u, bc)

    # Save solution to VTK
    res_file << u

    # Update previous solution
    u_n.assign(u)
