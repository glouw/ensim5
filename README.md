## ENSIM5

![](img/raylib2.png)

Ensim5 is a single threaded, branch free, allocation free, real time internal combustion engine audio DSP with 250Hz input/output control.

[https://www.youtube.com/watch?v=za0bQ6HqPRg](https://www.youtube.com/watch?v=za0bQ6HqPRg)

Ensim5's DPS engine simulates a thermofluidic otto-cycle with zero-dimensional isentropic flow
and one-dimensional pipe computational fluid dynamics.

## Build

Ensure raylib is installed (`pacman -S raylib`) and a c++20 compliant compiler like clang, then:

`make && ./gui`

## Controls

`1,2,3,4`: Throttle

`0`: Disables ignition

`Q,E`: Camera pan

`W,A,S,D`: Chamber select

`<,>`: Gear decrement/increment

### ATTRIBUTION

Ange Yaghi for the inspiration.
