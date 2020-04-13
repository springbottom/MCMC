
#include <iostream>
#include <string>
#include <fstream>
#include <vector>
#include <math.h>
#include <time.h>
#include <stdlib.h>

struct fdg_settings{
  int n;
  double epsilon;
  double A;
};

double torus_distance(double x1,double y1, double x2, double y2);

class fdg_loop : public MCMC_state{
private:
  double* state;
  double* energy; //represents all the bonds connected to that vertex...
  double* delta_energy;
  double* distances;
  fdg_settings* settings;
  int index;
  double new_x, new_y;
  double old_x, old_y;
  double* new_distances;
  double old_energy, new_energy;
  double total_delta_energy;

public:
  int num;
  fdg_loop(fdg_settings*);
  ~fdg_loop();
  int get_index(int i,int j);
  double update(double beta);
};
