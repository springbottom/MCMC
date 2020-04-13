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

int sweep(MCMC_settings* settings,
          MCMC_state* state,
          MCMC_stats* stats)
{
  double temp;
  for (int i = 0; i < (*settings).temp_N; i++){
    temp = (*settings).temp_min + i*((*settings).temp_max-(*settings).temp_min)/(double)((*settings).temp_N-1);
    (*stats).reset();
    for (int j = 0; j < (*settings).pre_iters; j++){
      (*state).update(1/temp);
    }
    for (int j = 0; j < (*settings).iters; j++){
      (*state).update(1/temp);
      (*stats).update(state,j+1);
    }
    (*stats).final(state,settings, 1/temp);
    //file << temp << "," << (*stats).final(state, 1/temp,(*settings).iters) << "\n";
  }

  //file.close();
  printf("Everything went well!\n");
  return 0;
}

MCMC_settings :: MCMC_settings(){

}

MCMC_settings :: MCMC_settings(int pre_iters, int iters, double temp_min, double temp_max, int temp_N, std::string file_name){
  (*this).pre_iters = pre_iters;
  (*this).iters = iters;
  (*this).temp_min = temp_min;
  (*this).temp_max = temp_max;
  (*this).temp_N = temp_N; //Actually we can just store like... an "iterator" I believe? representing the values of temp over this fopr loop.
  (*this).file_name = file_name;

  (*this).file = std::ofstream ((*this).file_name);
}

MCMC_settings :: ~MCMC_settings(){
  (*this).file.close();
}



//Welford stuff
void welford_state :: reset()
{
  old_x_n = 0;
  x_n = 0;
  M_2_n = 0;
}

void welford_state :: update(MCMC_state* state, int n)
{
  //Online algorithm for computing variance
  old_x_n = x_n;
  x_n   = old_x_n + ((*state).total_delta_energy - old_x_n)/(double)n;
  M_2_n = M_2_n + ((*state).total_delta_energy - old_x_n)*((*state).total_delta_energy - x_n);
}

void welford_state :: final(MCMC_state* state, MCMC_settings* settings, double beta)
{
  (*settings).file << beta << "," << beta*beta*M_2_n/(double)((*settings).iters * (*state).num) << "\n";
}
