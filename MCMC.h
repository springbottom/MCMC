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

template <class T_state, class T_stats> //Generic container for states?
int sweep(MCMC_settings* settings,
          T_state* init(void),
          T_stats* stats);

class welford_state
{
private:
  double old_x_n;
  double x_n;
  double M_2_n;
public:
  void reset();
  void update(double new_value, int n);
  double final(double beta);
};
