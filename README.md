## ENSIM5

![](img/raylib2.png)

Ensim5 is a single threaded, branch free, allocation free, real time internal combustion engine audio DSP with 250Hz input/output control.

[https://www.youtube.com/watch?v=za0bQ6HqPRg](https://www.youtube.com/watch?v=za0bQ6HqPRg)

## Build

Ensure raylib is installed (`pacman -S raylib`) and a c++20 compliant compiler like clang, then:

`make`
`./gui`

## Controls

1,2,3,4: Throttle
0: Disables ignition
Q,E: Camera pan
W,A,S,D: Chamber select
<,>: Gear decrement/increment

### Notes

Ensim5 improves on Ensim4 by exploring SIMD and cache locality for piston kinematics, isentropic flow, and computational fluid dynamics.

### ATTRIBUTION

Ange Yaghi for the inspiration.
