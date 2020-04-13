#include <string>
#include <fstream>

double rand_double();
int positive_modulo(int i, int n);

//Simple struct to hold my settings?
struct MCMC_settings
{
    int n;
    int pre_iters;
    int iters;
    double temp_min;
    double temp_max;
    int temp_N;
    std::string file_name;
};

//template <class T_state, class T_stats> //Generic container for states?
//int sweep(MCMC_settings*,T_state*,T_stats*);

template <class T_state, class T_stats> //Generic container for states?
int sweep(MCMC_settings* settings,
          T_state* state,
          T_stats* stats)
{
  std::ofstream file ((*settings).file_name);

  double temp;
  for (int i = 0; i < (*settings).temp_N; i++){
    temp = (*settings).temp_min + i*((*settings).temp_max-(*settings).temp_min)/(double)((*settings).temp_N-1);
    (*stats).reset();
    for (int j = 0; j < (*settings).pre_iters; j++){
      (*state).update(1/temp);
    }
    for (int j = 0; j < (*settings).iters; j++){
      (*stats).update((*state).update(1/temp),j+1);
    }
    file << temp << "," << (*stats).final(1/temp,(*state).num,(*settings).iters) << "\n";
  }

  file.close();
  printf("Everything went well!\n");
  return 0;
}




class welford_state
{
private:
  double old_x_n;
  double x_n;
  double M_2_n;
public:
  void reset();
  void update(double new_value, int n);
  double final(double beta, int num_samples, int num_particles);
};
