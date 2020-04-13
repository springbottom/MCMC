
#include <iostream>
#include <string>
#include <fstream>
#include <vector>
#include <math.h>
#include <stdlib.h>

struct fdg_settings{
  int n;
  double epsilon;
  double A;
};

double torus_distance(double x1,double y1, double x2, double y2);
