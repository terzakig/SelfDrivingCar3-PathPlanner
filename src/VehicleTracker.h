#ifndef VEHICLE_TRACKER_H
#define VEHICLE_TRACKER_H

struct VehicleTracker 
{
  int id;
  double x, y, s, d, vx, vy;
  double speed; // this should speed in the s presumably...
};

#endif