# cybertraining
Team 1 Project of the CyberTraining program at UMBC in 2018 (http://cybertraining.umbc.edu/)<br>
<br><b>Title:</b> Numerical Methods for Parallel Simulation of Diffusive Pollutant Transport from a Point Source <br>
<br><b>Team members:</b> Noah Sienkiewicz, Arjun Pandya, Tim Brown<br>
<br><b>Mentors:</b> Dr.Matthias K. Gobbert<br>
<br><b>Abstract</b>
<br>
In an interdisciplinary project combining Atmospheric Physics, High Performance Computing, and Big Data, we explore a numerical method for solving a physical system modeled by a partial differential equation. The application problem models the spread of pollution by a reaction-diffusion equation solved by the finite volume method. The numerical method is derived and tested on a known test problem in Matlab and then parallelized by MPI in C. We explore both closed and open systems of pollution, and show that the finite volume method is both mass conservative and has the ability to handle a point source modeled by the Dirac delta distribution. A parallel performance study confirms the scalability of the implementation to several compute nodes.
