#include <fstream>
#include <math.h>
#include <uWS/uWS.h>
#include <chrono>
#include <iostream>
#include <thread>
#include <vector>
#include "Eigen-3.3/Eigen/Core"
#include "Eigen-3.3/Eigen/QR"
#include "json.hpp"

#include "spline.h"
#include "VehicleTracker.h"

using namespace std;

const double LANE_SIZE = 4.0;         // in METERS!
const double TRACK_LENGTH = 6945.554; // in METERS!


// for convenience
using json = nlohmann::json;

// For converting back and forth between radians and degrees.
constexpr double pi() { return M_PI; }
double deg2rad(double x) { return x * pi() / 180; }
double rad2deg(double x) { return x * 180 / pi(); }

// Checks if the SocketIO event has JSON data.
// If there is data the JSON object in string format will be returned,
// else the empty string "" will be returned.
string hasData(string s) {
  auto found_null = s.find("null");
  auto b1 = s.find_first_of("[");
  auto b2 = s.find_first_of("}");
  if (found_null != string::npos) {
    return "";
  } else if (b1 != string::npos && b2 != string::npos) {
    return s.substr(b1, b2 - b1 + 2);
  }
  return "";
}

double distance(double x1, double y1, double x2, double y2)
{
	return sqrt((x2-x1)*(x2-x1)+(y2-y1)*(y2-y1));
}
int ClosestWaypoint(double x, double y, vector<double> maps_x, vector<double> maps_y)
{

	double closestLen = 100000; //large number
	int closestWaypoint = 0;

	for(int i = 0; i < maps_x.size(); i++)
	{
		double map_x = maps_x[i];
		double map_y = maps_y[i];
		double dist = distance(x,y,map_x,map_y);
		if(dist < closestLen)
		{
			closestLen = dist;
			closestWaypoint = i;
		}

	}

	return closestWaypoint;

}

int NextWaypoint(double x, double y, double theta, vector<double> maps_x, vector<double> maps_y)
{

	int closestWaypoint = ClosestWaypoint(x,y,maps_x,maps_y);

	double map_x = maps_x[closestWaypoint];
	double map_y = maps_y[closestWaypoint];

	double heading = atan2( (map_y-y),(map_x-x) );

	double angle = abs(theta-heading);

	if(angle > pi()/4)
	{
		closestWaypoint++;
	}

	return closestWaypoint;

}

// Transform from Cartesian x,y coordinates to Frenet s,d coordinates
vector<double> getFrenet(double x, double y, double theta, vector<double> maps_x, vector<double> maps_y)
{
	int next_wp = NextWaypoint(x,y, theta, maps_x,maps_y);

	int prev_wp;
	prev_wp = next_wp-1;
	if(next_wp == 0)
	{
		prev_wp  = maps_x.size()-1;
	}

	double n_x = maps_x[next_wp]-maps_x[prev_wp];
	double n_y = maps_y[next_wp]-maps_y[prev_wp];
	double x_x = x - maps_x[prev_wp];
	double x_y = y - maps_y[prev_wp];

	// find the projection of x onto n
	double proj_norm = (x_x*n_x+x_y*n_y)/(n_x*n_x+n_y*n_y);
	double proj_x = proj_norm*n_x;
	double proj_y = proj_norm*n_y;

	double frenet_d = distance(x_x,x_y,proj_x,proj_y);

	//see if d value is positive or negative by comparing it to a center point

	double center_x = 1000-maps_x[prev_wp];
	double center_y = 2000-maps_y[prev_wp];
	double centerToPos = distance(center_x,center_y,x_x,x_y);
	double centerToRef = distance(center_x,center_y,proj_x,proj_y);

	if(centerToPos <= centerToRef)
	{
		frenet_d *= -1;
	}

	// calculate s value
	double frenet_s = 0;
	for(int i = 0; i < prev_wp; i++)
	{
		frenet_s += distance(maps_x[i],maps_y[i],maps_x[i+1],maps_y[i+1]);
	}

	frenet_s += distance(0,0,proj_x,proj_y);

	return {frenet_s,frenet_d};

}

// Transform from Frenet s,d coordinates to Cartesian x,y
vector<double> getXY(double s, double d, vector<double> maps_s, vector<double> maps_x, vector<double> maps_y)
{
	int prev_wp = -1;

	while(s > maps_s[prev_wp+1] && (prev_wp < (int)(maps_s.size()-1) ))
	{
		prev_wp++;
	}

	int wp2 = (prev_wp+1)%maps_x.size();

	double heading = atan2((maps_y[wp2]-maps_y[prev_wp]),(maps_x[wp2]-maps_x[prev_wp]));
	// the x,y,s along the segment
	double seg_s = (s-maps_s[prev_wp]);

	double seg_x = maps_x[prev_wp]+seg_s*cos(heading);
	double seg_y = maps_y[prev_wp]+seg_s*sin(heading);

	double perp_heading = heading-pi()/2;

	double x = seg_x + d*cos(perp_heading);
	double y = seg_y + d*sin(perp_heading);

	return {x,y};

}

int getLane(double d)
{
  return (int)floor(d / LANE_SIZE);
}


// track distance. Now this is useful. It returns the distance between s2 and s1 in the track
// as "s2-s1"
double trackDistance(double s1, double s2)
{
  
  if (s2 > s1)
  {
    if (TRACK_LENGTH - s2 + s1 < s2 - s1) 
    {
      return TRACK_LENGTH - s2 + s1;
    }
    else
    {
      return s2 - s1;
    }
  }
  else
  {
    if (TRACK_LENGTH - s1 + s2 < s1 - s2) 
    {
      return -(TRACK_LENGTH - s1 + s2);
    }
    else
    {
      return -(s1 - s2);
    }
  }
  
  assert(false && "Unreachable region");
  
  return -999; // unreachable
}




// The planner
struct Planner
{
  // cost function constants
  const double COST_COLLISION = 1000000;
  const double COST_OFF_MIDDLE_LANE = 2;
  const double COST_JERK = 32;
  const double COST_MAX_ACCELERATION_VIOLATION = 28;
  const double COST_DEVIATION_FROM_SPEED_LIMIT = 3;
  const double COST_LANE_CHANGING = 1;
  
  // state constants
  /*static const int STATE_IDLE = 0; // car is idle (only at the beginning.
  static const int STATE_CHANGE_LANE_PREPARATION = 1; // preparing to change lane (possibly slow down to avoid jerk in the turn)
  static const int STATE_CHANGING_LANE = 2;           // changing lane
  static const int STATE_KEEP_LANE = 3;              // keep lane. 
  */
  // the current lane
  int current_lane;
  // the next lane
  int next_lane;
  
  // a copunter fopr low speed levels
  int low_speed_counter;
  
  // the speed limit in lane 0 (m/s)
  const double SPEED_LIMIT0 = 23.0; // go for it!
  // The speed limit in lane 1 (m/s)
  const double SPEED_LIMIT1 = 21.5;
  // the speed limit in lane 2 (m/s)
  const double SPEED_LIMIT2 = 19.0;
  
  // Maximum acceleration
  const double MAX_ACCELERATION = 10;
  // collision rtange
  const double COLLISION_RANGE = 3.0; 
  
  // Possible s-distances
  vector<double> s_distances;
  // possible end-speeds (in the s-axis)
  vector<double> s_end_speeds;
  
  // The top speed for the current lane
  double top_lane_speed; // speed limit in lane #1
  
  Planner()
  {
    s_distances = {20, 25, 30, 35, 50 , 55, 60, 70, 90 , 100, 110, 120/*, 160*/};
    
    current_lane = 1;
    next_lane = 1;
    top_lane_speed = SPEED_LIMIT1;
    s_end_speeds = {top_lane_speed * 1.2,
		    top_lane_speed * 1.1,
		    top_lane_speed,
		    top_lane_speed * 0.9, 
		    top_lane_speed * 0.8, 
		    top_lane_speed * 0.6, 
		    top_lane_speed * 0.5, 
		    top_lane_speed * 0.4, 
		    top_lane_speed * 0.3, 
		    top_lane_speed * 0.2,
		    top_lane_speed * 0.1,
		    top_lane_speed * 0.07
		    };
  
    low_speed_counter = 0;
  }
  
  
  pair<vector<double>, vector<double>> OptimalPlan1(const vector<VehicleTracker> &trackers,
						    const vector<double> &prev_x,
						    const vector<double> &prev_y,
						    double car_x,
						    double car_y,
						    double car_speed,
						    double car_yaw,
						    double car_s,
						    double car_d,
						    const vector<double> &map_s,   // the waypoint s-coordinates
						    const vector<double> &map_x,   // the waypoint s-coordinates
						    const vector<double> &map_y,   // the waypoint s-coordinates
						    double Dt = 0.02,       // Time step (20 ms)
						    int numSparse_sd_Points = 3 // number of sparse s-d points
						   )
  {
    // update the current lane
    current_lane = getLane(car_d);
    // get the top speed in this lane
    switch(current_lane)
    {
      case 0 :  
	top_lane_speed = SPEED_LIMIT0;
      break;
      case 1: 
	top_lane_speed = SPEED_LIMIT1;
      break;
      case 2: 
	top_lane_speed = SPEED_LIMIT2;
      break;
    }
    
    // update the low-speed counter
    //if (car_speed < 0.75 * top_lane_speed) 
    if (car_speed < 0.7 * top_lane_speed) 
    {
      low_speed_counter++;
    }
    else 
    {
      low_speed_counter = 0;
    }
    cout <<"Low speed counter : "<<low_speed_counter<<endl;
    
    
    
	
    
    
    double best_s_distance = -1;
    double bestT = -1;
    double best_end_speed = -1;
    double best_acc_s = -1;
    double minCost = 99999999999.9;
    
    	 
    // here we need to gradually increase velocity, so we seek a maximum s-distance in which we can veeeeeery gradually grow our speed
    // and at the same time avoid having somebody bumping on us from behind.The idea is to try to reach maximum speed, but without
    // accelerating fast.
    for (int i = 0; i < s_distances.size(); i++)
    {
      double s_distance = s_distances[i];
      //cout <<"s_distance : "<<s_distance<<endl;
      // future time (forward) for collision prediction
      double T_forward = 40 * Dt;
      double car_s_forward = car_s + s_distance;
      if (car_s_forward> TRACK_LENGTH) car_s_forward -= TRACK_LENGTH;
      
     s_end_speeds = {top_lane_speed * 1.2,
		    top_lane_speed * 1.1,
		    top_lane_speed,
		    top_lane_speed * 0.9, 
		    top_lane_speed * 0.8, 
		    top_lane_speed * 0.6, 
		    top_lane_speed * 0.5, 
		    top_lane_speed * 0.4, 
		    top_lane_speed * 0.3, 
		    top_lane_speed * 0.2,
		    top_lane_speed * 0.1,
		    top_lane_speed * 0.07
		    };
      // now check possible end-speeds in s
      for (int j = 0; j < s_end_speeds.size(); j++)
      {
	double end_speed = s_end_speeds[j];
	// cout <<"end_speed : "<<end_speed<<endl;
	// computing the time duration:
	double T = 2 * s_distance / (car_speed + end_speed);
	// Now computing the acceleration in s
	double acc_s = 2 * (s_distance - car_speed * T) / (T * T);
	    
	bool collision = false;
	double collision_distance = 0;
	// Checking for collisions
	for (int k = 0; k < trackers.size() && !collision; k++)
	{
	  VehicleTracker tracker = trackers[k];
	  if (tracker.speed == 0) continue;
	      
	  //cout << "Tracker's d "<<tracker.d<<" and lane : "<< getLane(tracker.d)<<endl;
	  
	  
	  if ( fabs(tracker.d - car_d) > 3 ) continue;
	  // check for immediatre collisioon before before anything...
	  if (fabs( trackDistance(tracker.s + tracker.speed * 5 * Dt, car_s + car_speed * 5 * Dt) ) < 8) // a few meters margin on the s-axis in current position
	  {
	    collision = true;
	  }
	  if ( getLane(tracker.d) != current_lane) continue;
	  // now project the s-location of the vehicle for time t = T
	  double tracker_s_forward = tracker.s + T_forward * tracker.speed;
	  if (tracker_s_forward > TRACK_LENGTH) tracker_s_forward -= TRACK_LENGTH;
	  if (trackDistance(tracker.s, car_s) < 0 && trackDistance(tracker_s_forward, car_s_forward) > 0) // the vehicle hit us from behind!
	  {
	    collision = true;  
	  }
	  if (trackDistance(tracker.s, car_s) > 0 && trackDistance(tracker_s_forward, car_s_forward) < 0) // we bumbed onto the vehicle!
	  {
	    collision = true;
	  }
	      
	 
	  collision_distance = fabs(trackDistance(tracker.s, car_s));
	}
	double cost = 0;    
	// penalize acceleration in s-axis only when keeping the same lane (i.e. trying to do lane change fast)
	if (current_lane == next_lane)  cost+= 4.5*abs(acc_s); // could punish shorter plans for potentially less jerk and acceleration of the spline
	else cost += 0.5 * fabs(acc_s) + T; // penalize the time duration of the manoeuvre
	if (top_lane_speed > end_speed) cost += fabs(top_lane_speed - end_speed);
	if (collision == true) cost += COST_COLLISION * collision_distance;
	if (cost < minCost)
	{
	  minCost = cost;
	  bestT = T;
	  best_acc_s = acc_s;
	  best_end_speed = end_speed;
	  best_s_distance = s_distance;
	      
	}
      }
	    
    }
	 
    // now switching state if necessary.
    // The car does not accelerate any more, although we still can go faster, so we maybe want to change lane...
    cout <<"Acceleration in s : "<<best_acc_s<<" and our next speed is : " << best_end_speed<<endl;
    //if (fabs(best_acc_s) < 0.1 && best_end_speed < 0.8 * SPEED_LIMIT && current_lane == next_lane)
    if (low_speed_counter > 135)
    {
      low_speed_counter = 0;
      cout << " We have been cruising in lower speed for too long. I think I want to CHANGE LANES... " <<endl;
      // Checking neighboring lanes 
      bool lane_found = false;
      int candidate_lane = current_lane;
      for (int lane_ofs = - 1; lane_ofs < 2 && !lane_found; lane_ofs++)
      {
	candidate_lane = current_lane + lane_ofs;
	cout <<"Candidate lane : "<<candidate_lane<<endl;
	// skip this step if the lane is invalid
	if (candidate_lane < 0 || candidate_lane > 2 || lane_ofs == 0) continue;
	// project the car's position in time (roughly...)
	double T_forward = 80 * Dt;
	double car_s_forward = car_s + car_speed * cos(pi() / 4) * T_forward;
	double car_d_forward = car_d + lane_ofs * car_speed * sin(pi() / 4) * T_forward;
	
	// check if the lane is clear
	bool lane_clear = true;
	for (int k = 0; k < trackers.size() && lane_clear; k++)
	{
	  VehicleTracker tracker = trackers[k];
	  // check for collision in d first (so that we can discard the ones thatr are far)
	  if ( fabs( tracker.d - (candidate_lane + 0.5)* LANE_SIZE  ) > 3 ) continue;
	  // now requiring that the track-distance of the vehicle from the car is at least than 40 meters
	  if ( fabs( trackDistance(car_s, tracker.s) ) < 8 ) lane_clear = false; // immediate collision
	  //project in time
	  double tracker_s_forward = tracker.s + T_forward * tracker.speed;
	  if (tracker_s_forward > TRACK_LENGTH) tracker_s_forward -= TRACK_LENGTH;
	  if (trackDistance(tracker.s, car_s) < 0 && trackDistance(tracker_s_forward, car_s_forward) > 0 && fabs(tracker.d - car_d_forward) < 3) // the vehicle hit us from behind!
	  {
	    lane_clear = false;  
	  }
	  if (trackDistance(tracker.s, car_s) > 0 && trackDistance(tracker_s_forward, car_s_forward) < 0 && fabs(tracker.d - car_d_forward) < 3) // we bumbed onto the vehicle!
	  {
	    lane_clear = true;
	  }
	      
	  
	}
	if (lane_clear) lane_found = true;
      }
      // do the lane change
      if (lane_found) next_lane = candidate_lane; 
      else next_lane = current_lane;
      cout <<"Changing lane from " <<current_lane << " to "<<next_lane<<endl;
    }
    
    cout <<"GREAT!!!!!!!!!!! OBTAINED AN OPTIMAL PLAN!"<<endl;
    cout <<"The min cost : "<<minCost<<endl;
	      cout <<"T : "<<bestT<<endl;
	      cout << "acceleration : "<<best_acc_s<<endl;
	      cout << "end speed : "<<best_end_speed<<endl;
	      cout <<"best distance : "<<best_s_distance<<endl;
    // Fit the x-y splines now...
    pair<tk::spline, tk::spline> splines;
    if (prev_x.size() < 2)
    {
      splines = GenerateXYSplines(car_x, 
					car_y,
					car_yaw,
					car_speed,
					car_s,
					car_d,
					best_s_distance,     // The total s-distance to travel
					bestT,     // The duration of the manoeuvre
					next_lane,      // The finishing lane
					map_s, // the waypoint s-coordinates
					map_x, // the waypoint s-coordinates
					map_y,  // the waypoint s-coordinates
					Dt,      // Time step (20 ms)
					numSparse_sd_Points // number of sparse s-d points
					);
    }
    else
    {
      splines = GenerateXYSplines(prev_x,
					prev_y,
					car_s,
					car_d,
					best_s_distance,     // The total s-distance to travel
					bestT,     // The duration of the manoeuvre
					next_lane,      // The finishing lane
					map_s, // the waypoint s-coordinates
					map_x, // the waypoint s-coordinates
					map_y,  // the waypoint s-coordinates
					Dt,        // Time step (20 ms)
					numSparse_sd_Points  // number of sparse s-d points
					);
    }
    
      // FINALLY, generate the trajectory points
    // 4. Generate enough trajectory points from 0 to T+13 * Dt
    vector<double> traj_x, traj_y;
  
    cout << "Best T : "<<bestT<<endl;
    
    for (int i = 0; Dt * i <= bestT + 13 * Dt; i++)
    {
      double t = i * Dt;
      double x_pos = splines.first(t);
      double y_pos = splines.second(t);
      traj_x.push_back( x_pos );
      traj_y.push_back( y_pos );
      //cout <<" Adding the "<<i<<"th"<<"point for t = "<<t<<endl;
    
      // Going through all tracked vehicles to detect collisions (I hope so at least....)
      for (int l = 0; l < trackers.size(); l++)
      {
	VehicleTracker tracker = trackers[l];
	     
	int vehicle_lane = (int)floor(tracker.d / LANE_SIZE);
	      
	if (vehicle_lane != (int)floor(car_d / LANE_SIZE) ) continue;
	
	double vehicle_s = tracker.s + t * tracker.speed;
	if (vehicle_s > TRACK_LENGTH) vehicle_s -= TRACK_LENGTH;
	vector<double> vcoords = getXY(vehicle_s, tracker.d, map_s, map_x, map_y);     
	      
	// Now if the distance is within collision range, increase the score by a collision penalty
	if ( distance(vcoords[0], vcoords[1], x_pos, y_pos) < COLLISION_RANGE ) 
	{
	  cout <<"Potential Collision with vehicle : "<<tracker.id<<endl;
	  
	}
      }
      
    }
   
     return pair<vector<double>, vector<double>>(traj_x, traj_y);
  }
  
  
  /*
  * Here's the function that generates the splines on the x and y axes.
  * The idea is to have a starting and end-position and as well as a time duration of the manoeuvre
  * and based on this input, generate a sparse sequence of s-d points, which can then be interpolated 
  * in the world frame (if we use the s-d coordinates, the transformed trajectory will end-up being discontinuous due to the linearization)
  * NOTE: The spline is parametrized by time.
  * NOTE: This is the version that generates uses the unfinished points from the previous plan execution by the simulator.
  */
  pair<tk::spline, tk::spline> GenerateXYSplines(const vector<double> &previous_x,
						 const vector<double> &previous_y,
						 double car_s,
						 double car_d,
						 double S,     // The total s-distance to travel
						 double T,     // The duration of the manoeuvre
						 int lane2,      // The finishing lane
						 const vector<double> &map_s, // the waypoint s-coordinates
						 const vector<double> &map_x, // the waypoint s-coordinates
						 const vector<double> &map_y,  // the waypoint s-coordinates
						 double Dt = 0.02,       // Time step (20 ms)
						 int numSparse_sd_Points = 3 // number of sparse s-d points
						 )
  {
    // NOTE !!!!!!! DONT FORGET TO GET ONLY A HANDFUL OF THE ORIGINAL PLAN'S POINTS!!!!!
    //              OTHERWISE WE WILL BE HAVING INXPLICABLE SIDE-EFFECTS ON-TRACK......
    vector<double> prev_x({previous_x[0], previous_x[1]});
    vector<double> prev_y({previous_y[0], previous_y[1]});
    
    
    // The points to be interpolated
    vector<double> points_x;
    vector<double> points_y;
    vector<double> points_t;
  
  
    // push these previous points back in the point-list
    double t_offset;
    for (int i = 0; i < prev_x.size(); i++)
    {
      t_offset = i * Dt;
      points_x.push_back(prev_x[i]);
      points_y.push_back(prev_y[i]);
      points_t.push_back(t_offset);
    }  
  
    // let's first estimate velocity at the end of the previous plan
    // NOTE: Assuming it's all concentrated in the s-axis
    double x_0 = prev_x[prev_x.size()-1], y_0 = prev_y[prev_y.size()-1];
    double x_1 = prev_x[prev_x.size()-2], y_1 = prev_y[prev_y.size()-2];
  
    double v1_s = sqrt( (x_0 - x_1) * (x_0 - x_1) + (y_0 - y_1) * (y_0 - y_1) ) / Dt;
  
    // 2. Now we need to generate 3 points in the s-axis with constant acceleration
  
    // Compute an acceleration in s (assuming that car_speed is mostly projected in the s-axis):
    double acc_s = 2 * (S - v1_s * T) / (T * T);
  
    // start and finish d-coordinates
    // NOTE: We don't reeeealy push the starting s-d coordinates, 'cause this will cause discontinuities.
    //       They are just used to generate the sparse intermediate s-d points.
    double s_start = car_s;
    double d_start = car_d,
	   d_finish = lane2 * LANE_SIZE  + LANE_SIZE * 0.5;
  
    // now computing the starting and finishing d-coordinates
    double s_inc = S / numSparse_sd_Points;
    double d_inc = (d_finish - d_start) / numSparse_sd_Points;
  
    double wp_s, wp_d, wp_x, wp_y, t; // cahce variables
    for (int i = 1; i <= numSparse_sd_Points; i++)
    {
    
      // generating the 3 remaining map_waypoints 
      wp_s = i * s_inc;           // relative s-0coordinates 
    
      wp_d = d_start + i * d_inc; // absolute d-coordinates
    
      // the respective time. We need to solve a quadratic...........
      // NOTE: ALL BEING well, we should always have a positive solution
      double a = 0.5 * acc_s, b = v1_s, c = -wp_s;
      double D = b*b - 4 * a * c;
      double t1 = (-b + sqrt(D)) / (2 * a),
	     t2 = (-b - sqrt(D)) / (2 * a);
    
      t = t1 > 0 ? t1 : t2;
      if (t == 0) t += Dt;
      // Now convert into world coordinates
      wp_s += s_start;
      if (wp_s > TRACK_LENGTH) wp_s -= TRACK_LENGTH;
    
      vector<double> wp_xy = getXY(wp_s, wp_d, map_s, map_x, map_y);
      wp_x = wp_xy[0],
      wp_y = wp_xy[1];
      
      // push everything in the interpolation lists
      points_x.push_back(wp_x);
      points_y.push_back(wp_y);
      points_t.push_back(t + t_offset);
    	   
    }
  
    // 3. FINALLY, we need to add another point to ensure that the car will be cruising parallel to the d-axis
    double v2_s = v1_s + acc_s * T;
    wp_s += 10 * v2_s * Dt;
  
    if (wp_s > TRACK_LENGTH) wp_s -= TRACK_LENGTH;
    
    vector<double> wp_xy = getXY(wp_s, wp_d, map_s, map_x, map_y);
    wp_x = wp_xy[0],
    wp_y = wp_xy[1];
    points_x.push_back(wp_x);
    points_y.push_back(wp_y);
    points_t.push_back(t + 10*Dt + t_offset);
  
    // All being well, we vcan fit two splines in the x and the y-axis respectively
    tk::spline spline_x;
    spline_x.set_points(points_t,points_x);    

    tk::spline spline_y;
    spline_y.set_points(points_t,points_y);    

    /*
    // 4. Generate enough trajectory points from 0 to T+13 * Dt
    vector<double> traj_x, traj_y;
  
    for (int i = 0; Dt * i <= T + ( 10 + previous_x.size() ) * Dt; i++)
    {
      t = i * Dt;
      traj_x.push_back( spline_x(t) );
      traj_y.push_back( spline_y(t) );
    
      //cout <<" Adding the "<<i<<"th"<<"point for t = "<<t<<endl;
    }
    */
    return pair<tk::spline, tk::spline>(spline_x, spline_y);
  }
  
  
  /*
  * Here's the function that generates the points that will be interpolated by a spline.
  * The idea is to have a starting and end-position and as well as starting and end-velocity 
  * and based on this input, generate a sparse sequence of s-d points, which can then be interpolated 
  * in the world frame (if we use the s-d coordinates, the transformed trajectory will end-up being discontinuous due to the linearization)
  * NOTE: The spline is parametrized by time.
  * NOTE: This is the version that generates two points (what a trick) from the car orientation in order to force tangent matching
  *       numerically (the spline tool does not do that explicitly)
  */
  pair<tk::spline, tk::spline> GenerateXYSplines(double car_x, 
						 double car_y,
						 double car_yaw,
						 double car_speed,
						 double car_s,
						 double car_d,
						 double S,     // The total s-distance to travel
						 double T,     // The duration of the manoeuvre
						 int lane2,      // The finishing lane
						 vector<double> map_s, // the waypoint s-coordinates
						 vector<double> map_x, // the waypoint s-coordinates
						 vector<double> map_y,  // the waypoint s-coordinates
						 double Dt = 0.02,       // Time step (20 ms)
						 int numSparse_sd_Points = 3 // number of sparse s-d points
						 )
  {
    // The points to be interpolated
    vector<double> points_x;
    vector<double> points_y;
    vector<double> points_t;
  
  
    // generating a direction vector 
    double u_x = cos(car_yaw), 
	   u_y = sin(car_yaw);
    // 1. Now get a couple of points in the inverse direction along u.
    /* NOTE: Let's assume that the car has been travelling at the announced speed for T = 0.04 seconds.
    *        
    */
  

  
   // push these points in the lists
    if (car_speed > 0) 
    {
      // Move twice back at the given velocity
      points_x.push_back(car_x - car_speed * 2 * Dt * u_x);
      points_y.push_back(car_y - car_speed * 2 * Dt * u_y);
      points_t.push_back(0);
    
      // move once back at given velocity
      points_x.push_back(car_x - car_speed * Dt * u_x);
      points_y.push_back(car_y - car_speed * Dt * u_y);
      points_t.push_back(Dt);
    
      // Add the current position now...
      points_x.push_back(car_x);
      points_y.push_back(car_y); 
      points_t.push_back(2 * Dt);
  
    }
    else
    {
      // In this case (zero car speed) we move slightly forward...
      points_x.push_back(car_x);
      points_y.push_back(car_y);
      points_t.push_back(0);
    
      // Move forward at a small speed
      points_x.push_back(car_x + 0.001 * Dt * u_x);
      points_y.push_back(car_y + 0.001 * Dt * u_y);
      points_t.push_back(Dt);
    
      // And again, move a bit forward at a small speed
      points_x.push_back(car_x + 0.001 * (2 * Dt) * u_x);
      points_y.push_back(car_y + 0.001 * (2 * Dt) * u_y);
      points_t.push_back(2 * Dt);
    
    }
  
  
    double t_offset = points_t[points_t.size()-1]; // remmebr, we pushed some points!!!!
    // 2. Now we need to generate 3 points in the s-axis with constant acceleration
  
    // Compute an acceleration in s (assuming that car_speed is mostly projected in the s-axis):
    double acc_s = 2 * (S - car_speed * T) / (T * T);
    cout <<"Estimated acceleration : "<<acc_s<<endl;
    cout <<"The respective duration : "<<T<<endl;
    // start and finish d-coordinates
    double s_start = car_s;
    double d_start = car_d,
	   d_finish = lane2 * LANE_SIZE  + LANE_SIZE * 0.5;
  
    // now computing the starting and finishing d-coordinates
    double s_inc = S / numSparse_sd_Points;
    double d_inc = (d_finish - d_start) / numSparse_sd_Points;
  
    double wp_s, wp_d, wp_x, wp_y, t; // cahce variables
    for (int i = 0; i < numSparse_sd_Points; i++)
    {
    
      // generating the 3 remaining map_waypoints 
      wp_s = i * s_inc;
      wp_d = d_start + i * d_inc;
    
      // the respective time. We need to solve a quadratic...........
      // NOTE: ALL BEING well, we should always have a positive solution
      double a = 0.5 * acc_s, b = car_speed, c = -wp_s;
      double D = b*b - 4 * a * c;
      double t1 = (-b + sqrt(D)) / (2 * a),
	     t2 = (-b - sqrt(D)) / (2 * a);
    
      t = t1 > 0 ? t1 : t2;
      if (t == 0) t += Dt;
    
      // Now convert into world coordinates
      wp_s += s_start;
      if (wp_s > TRACK_LENGTH) wp_s -= TRACK_LENGTH;
      vector<double> wp_xy = getXY(wp_s, wp_d, map_s, map_x, map_y);
      wp_x = wp_xy[0],
      wp_y = wp_xy[1];
      // push everything in the interpolation lists
      points_x.push_back(wp_x);
      points_y.push_back(wp_y);
      points_t.push_back(t + t_offset);
	   
    }
    // 3. FINALLY, we need to add another point to ensure that the car will be crusing paralklel to the d-axis
    double v2_s = car_speed + acc_s * T;
    wp_s += 10 * v2_s * Dt;
    if (wp_s > TRACK_LENGTH) wp_s -= TRACK_LENGTH;
    
    vector<double> wp_xy = getXY(wp_s, wp_d, map_s, map_x, map_y);
    wp_x = wp_xy[0],
    wp_y = wp_xy[1];
    points_x.push_back(wp_x);
    points_y.push_back(wp_y);
    points_t.push_back(t + 10*Dt + t_offset);
  
    // All being well, we can fit trwo splines in the x and the y-axis respectively
    tk::spline spline_x;
    spline_x.set_points(points_t,points_x);    

    tk::spline spline_y;
    spline_y.set_points(points_t,points_y);    

    // 4. Generate enough trajectory points from 0 to T+13 * Dt
    /* vector<double> traj_x, traj_y;
  
    for (int i = 0; Dt * i <= T + 13 * Dt; i++)
    {
      t = i * Dt;
      traj_x.push_back( spline_x(t) );
      traj_y.push_back( spline_y(t) );
      //cout <<" Adding the "<<i<<"th"<<"point for t = "<<t<<endl;
    
    }
    */
    return pair<tk::spline, tk::spline>(spline_x, spline_y);
  
  }
  
  
  
  
  
};


int main() 
{
  uWS::Hub h;

  // Load up map values for waypoint's x,y,s and d normalized normal vectors
  vector<double> map_waypoints_x;
  vector<double> map_waypoints_y;
  vector<double> map_waypoints_s;
  vector<double> map_waypoints_dx;
  vector<double> map_waypoints_dy;

  // Waypoint map to read from
  string map_file_ = "../data/highway_map.csv";
  // The max s value before wrapping around the track back to 0
  double max_s = 6945.554;

  ifstream in_map_(map_file_.c_str(), ifstream::in);

  string line;
  while (getline(in_map_, line)) {
  	istringstream iss(line);
  	double x;
  	double y;
  	float s;
  	float d_x;
  	float d_y;
  	iss >> x;
  	iss >> y;
  	iss >> s;
  	iss >> d_x;
  	iss >> d_y;
  	map_waypoints_x.push_back(x);
  	map_waypoints_y.push_back(y);
  	map_waypoints_s.push_back(s);
  	map_waypoints_dx.push_back(d_x);
  	map_waypoints_dy.push_back(d_y);
  }

  
  Planner planner;
  vector<VehicleTracker> trackers;
  
  h.onMessage([&planner, &trackers, &map_waypoints_x,&map_waypoints_y,&map_waypoints_s,&map_waypoints_dx,&map_waypoints_dy](uWS::WebSocket<uWS::SERVER> ws, char *data, size_t length,
                     uWS::OpCode opCode) {
    // "42" at the start of the message means there's a websocket message event.
    // The 4 signifies a websocket message
    // The 2 signifies a websocket event
    //auto sdata = string(data).substr(0, length);
    //cout << sdata << endl;
    if (length && length > 2 && data[0] == '4' && data[1] == '2') {

      auto s = hasData(data);

      if (s != "") {
        auto j = json::parse(s);
        
        string event = j[0].get<string>();
        
        if (event == "telemetry") {
          // j[1] is the data JSON object
          
        	// Main car's localization Data
          	double car_x = j[1]["x"];
          	double car_y = j[1]["y"];
          	double car_s = j[1]["s"];
          	double car_d = j[1]["d"];
          	double car_yaw_deg = j[1]["yaw"];
		double car_yaw = deg2rad(car_yaw_deg);
          	double car_speed_MPH = j[1]["speed"];
		double car_speed = 0.44704 * car_speed_MPH;
	  
          	// Previous path data given to the Planner
          	auto previous_path_x = j[1]["previous_path_x"];
          	auto previous_path_y = j[1]["previous_path_y"];
          	//cout <<"Previous path data-x size : "<<previous_path_x<<endl;
		// Previous path's end s and d values 
          	double end_path_s = j[1]["end_path_s"];
          	double end_path_d = j[1]["end_path_d"];

		cout <<"Car s: "<<car_s<<endl;
		cout <<"Car d: "<<car_d<<endl;
		
          	// Sensor Fusion Data, a list of all other cars on the same side of the road.
          	auto sensor_fusion = j[1]["sensor_fusion"];
		auto map = j[1]["map"];
		//cout <<"Sensor fusion 0 : "<<sensor_fusion[1][0]<<endl;
		json msgJson;

          	vector<double> next_x_vals;
          	vector<double> next_y_vals;

		pair<vector<double>, vector<double>> plan;
		
		// update the trackers
		//cout <<"^&R*&R* Sensor fusion : "<<sensor_fusion<<endl;
		if (trackers.size() == 0)
		{
		  VehicleTracker tracker; // cache strcuture
		  for (int j = 0; j < sensor_fusion.size(); j++)
		  {
		    tracker.id = (int)sensor_fusion[j][0];
		    tracker.x = (double)sensor_fusion[j][1];
		    tracker.y = (double)sensor_fusion[j][2];
		    tracker.vx = (double)sensor_fusion[j][3];
		    tracker.vy = (double)sensor_fusion[j][4];
		    tracker.s = (double)sensor_fusion[j][5];
		    tracker.d = (double)sensor_fusion[j][6];
		    tracker.speed = sqrt(tracker.vx * tracker.vx + tracker.vy * tracker.vy);
		    
		    
		    trackers.push_back(tracker);
		    /* cout <<"Updated tracvker : "<<trackers[j].id<<endl;
		  cout <<"Tracker x : "<<trackers[j].x<<endl;
		  cout <<"Tracker y : "<<trackers[j].y<<endl;
		  cout <<"Tracker vx : "<<trackers[j].vx<<endl;
		  cout <<"Tracker vy : "<<trackers[j].vy<<endl;
		  cout <<"Tracker s : "<<trackers[j].s<<endl;
		  cout <<"Tracker d : "<<trackers[j].d<<endl;
		  */
		  }
		}
		else
		{
		  for (int j = 0; j < sensor_fusion.size(); j++)
		  {
		    //trackers[j].id = (int)sensor_fusion[j][0];
		    trackers[j].x = (double)sensor_fusion[j][1];
		    trackers[j].y = (double)sensor_fusion[j][2];
		    trackers[j].vx = (double)sensor_fusion[j][3];
		    trackers[j].vy = (double)sensor_fusion[j][4];
		    trackers[j].s = (double)sensor_fusion[j][5];
		    trackers[j].d = (double)sensor_fusion[j][6];
		    trackers[j].speed = sqrt(trackers[j].vx * trackers[j].vx + trackers[j].vy * trackers[j].vy);
		    
		   /*cout <<"Updated tracvker : "<<trackers[j].id<<endl;
		  cout <<"Tracker x : "<<trackers[j].x<<endl;
		  cout <<"Tracker y : "<<trackers[j].y<<endl;
		  cout <<"Tracker vx : "<<trackers[j].vx<<endl;
		  cout <<"Tracker vy : "<<trackers[j].vy<<endl;
		  cout <<"Tracker s : "<<trackers[j].s<<endl;
		  cout <<"Tracker d : "<<trackers[j].d<<endl;
		  */
		    
		  }
		  
		  
		} 
		
	
		plan = planner.OptimalPlan1(trackers,
						   previous_path_x,
						   previous_path_y,
						   car_x,
						   car_y,
						   car_speed,
						   car_yaw,
						   car_s,
						   car_d,
						   map_waypoints_s,   // the waypoint s-coordinates
						   map_waypoints_x,   // the waypoint s-coordinates
						   map_waypoints_y,   // the waypoint s-coordinates
						   0.02,       // Time step (20 ms)
						   3 // number of sparse s-d points
						  );
		
		cout <<"New plan length : "<< (plan.first).size() << endl;
		
		for (int i = 0; i< plan.first.size(); i++)
		{
		  next_x_vals.push_back(plan.first[i]);
		  next_y_vals.push_back(plan.second[i]);
		  
		}
		
          	msgJson["next_x"] = next_x_vals;
          	msgJson["next_y"] = next_y_vals;

          	auto msg = "42[\"control\","+ msgJson.dump()+"]";

          	//this_thread::sleep_for(chrono::milliseconds(1000));
          	ws.send(msg.data(), msg.length(), uWS::OpCode::TEXT);
          
        }
      } else {
        // Manual driving
        std::string msg = "42[\"manual\",{}]";
        ws.send(msg.data(), msg.length(), uWS::OpCode::TEXT);
      }
    }
  });

  // We don't need this since we're not using HTTP but if it's removed the
  // program
  // doesn't compile :-(
  h.onHttpRequest([](uWS::HttpResponse *res, uWS::HttpRequest req, char *data,
                     size_t, size_t) {
    const std::string s = "<h1>Hello world!</h1>";
    if (req.getUrl().valueLength == 1) {
      res->end(s.data(), s.length());
    } else {
      // i guess this should be done more gracefully?
      res->end(nullptr, 0);
    }
  });

  h.onConnection([&h](uWS::WebSocket<uWS::SERVER> ws, uWS::HttpRequest req) {
    std::cout << "Connected!!!" << std::endl;
  });

  h.onDisconnection([&h](uWS::WebSocket<uWS::SERVER> ws, int code,
                         char *message, size_t length) {
    ws.close();
    std::cout << "Disconnected" << std::endl;
  });

  int port = 4567;
  if (h.listen(port)) {
    std::cout << "Listening to port " << port << std::endl;
  } else {
    std::cerr << "Failed to listen to port" << std::endl;
    return -1;
  }
  h.run();
}
















































































