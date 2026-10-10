# Leader-Follower-STL

MATLAB code for planning the motion of a follower robot in a leader–follower problem with signal temporal logic (STL) specifications. The specifications are encoded as mixed-integer constraints and the whole trajectory of the follower is computed by solving a mixed-integer linear program.

**For an explanation of the method and the results, see the project page:**
**https://mhyoosefian.github.io/projects/leader-follower-stl/**

## Problem

Two agents take part in a surveillance mission. The leader starts at $(0, -2)$ and moves counterclockwise on a circle of radius 2 m centered at the origin, with an angular velocity of 10°/s. The follower starts at $(-1, -3)$ and must:

- always stay within 2 m of the leader, so that the team does not separate;
- never come closer than 0.5 m to the leader, to avoid a collision;
- eventually visit the goal region $[-4.5, -3] \times [-0.5, 0.5]$ within the first 60 seconds;
- move with velocities $|u_x|, |u_y| \le 0.4$ m/s.

Both agents are modeled as points, and the follower moves with $x_{t+1} = x_t + u_x\,\Delta t$ and $y_{t+1} = y_t + u_y\,\Delta t$.

## Formulation

The separation requirements are written with the infinity norm, which is conservative with respect to the Euclidean distance: $\|p_f - p_l\|_\infty \le \sqrt{2}$ guarantees a distance of at most 2 m, and $\|p_f - p_l\|_\infty \ge 0.5$ guarantees at least 0.5 m. The full specification is

```math
\varphi = F_{[0,60]}\big(p_f \in \text{goal}\big) \wedge G_{[0,T]}\big(\|p_f - p_l\|_\infty \le \sqrt{2}\big) \wedge G_{[0,T]}\big(\|p_f - p_l\|_\infty \ge 0.5\big),
```

where $F$ means *eventually* and $G$ means *always*. The first separation requirement is a conjunction of linear constraints and is added directly. The goal region and the collision avoidance (a disjunction: the follower must be to the left, to the right, above or below the leader) are encoded with binary variables and the big-M method, and the logical and temporal operators become linear inequalities between these binary variables. Among all trajectories that satisfy the specification, the code finds the one with the least control effort, $\sum_t \big(|u_x(t)| + |u_y(t)|\big)$.

## Requirements

- MATLAB
- [YALMIP](https://yalmip.github.io/)
- [MOSEK](https://www.mosek.com/) (any mixed-integer solver supported by YALMIP can be used by changing the `solver` option)

## How to use the code

Run `code.m`. It solves the problem and then shows:

1. an animation of the leader and the follower, with the region the follower must stay in, the keep-out zone and the goal region (set `saveVideo = 1` to also save it as an MP4 file);
2. the velocity commands of the follower;
3. the distance between the two agents, with the two limits.

The problem parameters (the leader's path, the initial positions, the goal region, the time horizon and the velocity bounds) are set at the top of `code.m`.

## Repository structure

```
code.m                     problem definition, optimization and plots
animateLeaderFollower.m    animation of the solution
```
