# Goo Milestone 1
Michael Brenan (mb45487) & John Herrick (jrh5735)

## Build Instructions

No change from the default code - use the provided build instructions. We recommend building with
release mode to make the simulation not slow down with moderate numbers of particles. Note that the
submitted code comes with all of the third party dependencies already included.

```
mkdir build
cd build
cmake -DCMAKE_BUILD_TYPE=Release ..
make -j `nproc`
```

## Creative Component / Extension

None for this milestone; we'll do more next milestone ;)

## Discussion: Floor Force

We implemented a simple, conservative floor force based off of the potential energy function

V(q) = -(k/2)(-0.5 - q.y)^2     (if q.y < -0.05)
     = 0                  (otherwise)

Where q is the position of a particle; k is a strength constant which affects how quickly the particle rebounds off of
the floor; we set it to 25000 by default (so quite strong). Differentiating this value function yields a force/Jacobian and
Hessian of

{F(q)}(dq) = k(-0.5 - q.y)
{dF(q)}(dq)(dq2) = -k
