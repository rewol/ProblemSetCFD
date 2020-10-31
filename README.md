REVISION DATE: 31.10.2020
Computational Fluid Dynamics Basics with Applications by John D. Anderson
Chapter 10 Problem
Supersonic flow over a flat plate.
Refer to the book mor the methodology.

Altitude class calculates freestream values up to 6000m
Field2d class creates a flow field and sets dx, dy using the input arguments.
Solver class creates a solution field with Vector U, E and F.

Place the .h and .cpp files in the same directory as main.cpp and run the code.

1. Only MACCORMACK' s method is implemented in the code.
2. Only constant temperature on the wall is implemented. It'd be a nice challenge to implement adiabatic wall.
Make sure to edit boundary conditions and vector calculations.
3. MACCORMACK method is explicit so timestep calculation is limited with CFL stability criteria.
4. u, v, rho, T and p are extracted to text files in a format that is required by gnuplot contour plot.
5. Other types of plots can be extracted from the solution field as you see fit.
6. There is no massflow check at the inlet and outlet. It should definetely be implemented.
7. It is definetely possible to optimize the code by performing highlevel move / copy techniques. In my computer it took 417 seconds to converge for a 121 by 121 grid.
After 121 nodes, solution becomes unstable and "explode". So for grid independence for this particular problem in the book, a more detail look into time step function is in order.

Play with the code as you want, use it as a basis for different problems and make experiments.
This is a very very very simple code and it has many many weaknesses.

Have fun and be happy.

Refik
