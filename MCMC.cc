/*
A general MCMC script.
*/

#include "MCMC.h"

//Random useful functions
double rand_double(){
  return rand()/(double)RAND_MAX;
}

int positive_modulo(int i, int n) {
    return (i % n + n) % n;
}

//Welford stuff
void welford_state :: reset()
{
  old_x_n = 0;
  x_n = 0;
  M_2_n = 0;
}

void welford_state :: update(double new_value, int n)
{
  //Online algorithm for computing variance
  old_x_n = x_n;
  x_n   = old_x_n + (new_value - old_x_n)/(double)n;
  M_2_n = M_2_n + (new_value - old_x_n)*(new_value - x_n);
}

double welford_state :: final(double beta, int num_samples, int num_particles)
{
  return beta*beta*M_2_n/(double)(num_samples * num_particles);
}
