/*
A general MCMC script.
*/
#include <string>


template <class T_state> //Generic container for states?
template <class T_stats> //Generic container for the stats module

//Simple struct to hold my settings?
struct MCMC_settings
{
    int n; int pre_iters; int iters;
    double temp_min; double temp_max; int temp_N;
    string file_name;
}

//The actual sweeping function
int sweep(MCMC_settings settings,
          T_state* init(void),
          double step(double,*T_state),
          ){

  ofstream file (file_name);
  T_state* state = init();

  double temp;
  for (int i = 0; i < N; i++){
    for (int j = 0; j < pre_iters; j++){
      step(1/temp,*state);
    }
    temp = temp_min + i*(temp_max-temp_min)/(double)(N-1);
    step(1/temp,*state);
  }



  my_file.close();


  printf("Everything went well!\n");
  return 0;
}

//Simple class to do the welford online variance calculation.
class welford_state
{
  double old_x_n;
  double x_n;
  double M_2_n;

  int update(double new_value, int n){
    //Online algorithm for computing variance
    this.old_x_n = x_n;
    this.x_n   = this.old_x_n + (new_value - this.old_x_n)/(double)n;
    this.M_2_n = this.M_2_n + (new_value - this.old_x_n)*(new_value - this.x_n);

    return 0;
  }

  double final(double beta){
    return beta*beta*this.M_2_n;
  }
}
