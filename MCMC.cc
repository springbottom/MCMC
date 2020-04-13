/*
A general MCMC script.
*/

#include "MCMC.h"

//The actual sweeping function
template <class T_state, class T_stats> //Generic container for states?
int sweep(MCMC_settings* settings,
          T_state* init(void),
          T_stats* stats)
{
  std::ofstream file ((*settings).file_name);
  T_state* state = init();

  double temp;
  for (int i = 0; i < (*settings).temp_N; i++){
    temp = (*settings).temp_min + i*((*settings).temp_max-(*settings).temp_min)/(double)((*settings).temp_N-1);
    (*stats).reset();
    for (int j = 0; j < (*settings).pre_iters; j++){
      (*state).step(1/temp);
    }
    for (int j = 0; j < (*settings).iters; j++){
      (*stats).update((*state).step(1/temp));
    }
    file << temp << "," << (*stats).final();
  }

  file.close();
  printf("Everything went well!\n");
  return 0;
}

//Simple class to do the welford online variance calculation.
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

double welford_state :: final(double beta)
{
  return beta*beta*M_2_n;
}

double rand_double(){
  return random()/(double)RAND_MAX;
}

inline int positive_modulo(int i, int n) {
    return (i % n + n) % n;
}
