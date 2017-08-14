# SelfDrivingCar3-PathPlanner
Real-time path planning for a car moving in a 3-lane track with traffic.

## Problem
This path planner navigates a car through traffic in a 3-lane track. The traffic simulator can be downloaded [here](https://github.com/udacity/self-driving-car-sim/releases/tag/T3_v1.2). The car should be able to avoid collisions with other cars and produce motion with minimal jerk, acceleration and speed below the speed limit. Furthermore, the car should perform lane changing quickly enough to avoid being classified as "off-lane". The objective is to drive for at least 4.32 miles without any of the aforementioned alerts showing-up on the simulator.  

## Overview
The overall principle of plannning here is to generate multiple different splines that reresent a lane change (or keep) and choose the one that minimizes the criteria stated above. The behavioral level is simple: It decides on whether the car should keep of change lane, based on the latest speed levels in the current lane. The trajectory design involves two splines, one for the x-coordinate and one for the y-coordinate of the planned trajectory in world coordinates by interpolating a sparse setr of s-d values, produced using a constant acceleration model for the s-axis motion of the car. More details in the [writeup pdf](https://github.com/terzakig/SelfDrivingCar3-PathPlanner/blob/master/writeup_report.pdf).

## Performance
This is a plannner which penalizes collisions, jerk, acceleration, deviation from speed limit in a cost function and therefore treats them as _soft_ constraints. In other words, with the exception of collisions, these quantities may sometimes "jump" above the given boundaries at any time during regular cruise, but that doesn't happen very often. However, in overall, the car achieves the 4.32 miles milestone rather easily in a generally safe (no collisions) manner. I would argue that the collision detection mechanism can be greatly improved, but for this application it generally does the job.   
