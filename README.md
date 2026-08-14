### Ensim5

1,2,3,4 for throttle control.

Q,E to cycle popups.

![]("img/img.png")

Ensim5 improves on Ensim4 by exploring SIMD and cache locality for piston kinematics and isentropic flow.

Intel(R) Core(TM) i7-8665U CPU @ 1.90GHz perf-stat.
48000 audio samples of a 36-chamber 4-piston engine and one dimensional CFD pipe generated in 0.2 seconds:

```
            0  context-switches:u       #  0.0    cs_per_second
            0  cpu-migrations:u         #  0.0    migrations_per_second
    2,376,758  L1-dcache-load-misses:u  #  0.7 %  l1d_miss_rate
       83,108  branch-misses:u          #  0.1 %  branch_miss_rate
1,627,262,496  instructions:u           #  2.2    insn_per_cycle
  362,993,860  dTLB-loads:u             #  0.0 %  dtlb_miss_rate
```

Unlike Ensim4, there is no bulk chamber resevoir to combine several SIMD lanes which does harm pipe audio fidelity.
