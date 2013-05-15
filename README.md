Only simulation.sce is working!
Simulation-parameters can be found at the beginning.
Exposition is defined at the bottom.
Execute with exec function in Scilab.

What the code dose:
It simlates the motion of charged point-particles at relativistic speeds (= near the speed of light),
using leap-frog-integration.

How?
First you need an exposition. I chose a disc-formed distribuition of random particles.
The accileration of this particles due to the forces they act on eachother is calculated.
At the beginning of each simulation step, it calculates the velocity + (accileration * timstep / 2), which will be the
velocity used for all further calcualtions during this step.
Now each particle doesn't see the particles as they are, but as they were since they sent out the
"light" that reaches them. (the field of charge travels with the speed of light.)
Therefore, we take a look back in our timeline and choose the particles which are in the right distance
so the particle can "feel" it or rather the forces of it's field.
Since we don't want to take any particle twice, we store all picked particles in a boolean-matrix
so we don't pick them again. (It has a false-diagonl, since the diagonal are the particles themselfs.)
Then the field each particle experiences by it's own particle-distribution is calculated. The force,
acctually combined of electrical, inductive and magnetic force, is computed through a neat formula
thanks to a paper on the topic (link in the Wiki), which only needs the electrical field of each particle.
All these force are again Lorentz-transformed to the system of the current particle.
(Before we transformed the particles acting on our current particle.)
The forces are then added up, so we get the force currently applied to the particle.
But since we are at relativistic speeds, we can't simply divide by the particles mass,
since it's mass will be much greater!
Now the "correction"-function calculates us the neccassery factor to get the correct acceleration.
The acceleration is added to the current velocity, but only for half a timestep (according to leap-frog-integration)!
This is the velocity we save, while we also save the acceleration, so we can calculate the speed for our next iteration.
