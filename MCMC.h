

#include <string>
#include <fstream>

double rand_double();
int positive_modulo(int i, int n);

//Simple struct to hold my settings?
class MCMC_settings
{
public:
  MCMC_settings();
  MCMC_settings(int pre_iters, int iters, double temp_min, double temp_max, int temp_N, std::string file_name);
  ~MCMC_settings();
  int pre_iters;
  int iters;
  double temp_min;
  double temp_max;
  int temp_N;
  std::string file_name;
  std::ofstream file;
};

//A state class?!
class MCMC_state
{
private:

public:
  int num; //num particles
  double total_delta_energy;
  void update(double beta){};
};

//A stats class
class MCMC_stats
{
private:

public:
  void reset(){};
  void update(MCMC_state*,int){};
  void final(MCMC_state*, MCMC_settings*, double beta){};
};


//template <class T_state, class T_stats> //Generic container for states?
int sweep(MCMC_settings* settings,
          MCMC_state* state,
          MCMC_stats* stats);


class welford_state : public MCMC_stats
{
private:
  double old_x_n;
  double x_n;
  double M_2_n;
public:
  void reset();
  void update(MCMC_state*, int n);
  void final(MCMC_state*, MCMC_settings*, double beta);
};
