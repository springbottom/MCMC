#include "MCMC.h"
#include "fdg.h"

int main(int argc, char *argv[]){
  time_t t;
  srand((unsigned) time(&t));

  MCMC_settings settings;
  fdg_settings settings2;

  if (argc == 10){
    settings2.n        = atoi(argv[1]);         //how many nodes in the thingo
    settings2.epsilon  = strtod(argv[2],NULL);  //the radius of hopping, probably keep it relatively small?
    settings2.A        = strtod(argv[3],NULL);  //H = \sum(bonds) |x_i - y_i| - A \sum(all) |x_i - y_i|

    settings.pre_iters = atoi(argv[4]);         //
    settings.iters     = atoi(argv[5]);         //
    settings.temp_min  = strtod(argv[6],NULL);  //
    settings.temp_max  = strtod(argv[7],NULL);  //
    settings.temp_N    = atoi(argv[8]);         //
    settings.file_name = argv[9];
  }
  else
  {
    printf("There are not enough arguments you silly!\n");
    settings2.n = 8;
    settings2.epsilon = 0.1;
    settings2.A = 0.5;

    settings.pre_iters = 10;
    settings.iters = 10;
    settings.temp_min = 0.2;
    settings.temp_max = 1.5;
    settings.temp_N = 10;
    settings.file_name = "test";
  }

  //std::cout << settings.file_name << std::endl;

  fdg_loop x = fdg_loop(&settings2);
  welford_state y = welford_state();
  sweep(&settings,&x,&y);

  return 0;
}

double torus_distance(double x1,double y1, double x2, double y2){
  //assuming |x_i|,|y_i|<1
  double delta_x = x2-x1;
  double delta_y = y2-y1;
  (delta_x > 0.5  ? delta_x-- : delta_x);
  (delta_x < -0.5 ? delta_x++ : delta_x);
  (delta_y > 0.5  ? delta_y-- : delta_y);
  (delta_y < -0.5 ? delta_y++ : delta_y);
  return pow(pow(delta_x,2) + pow(delta_y,2),0.5);
}

void torus_mod(double* x, double* y){
  if (*x < 0) (*x)++;
  if (*x > 1) (*x)--;
  if (*y < 0) (*y)++;
  if (*y > 1) (*y)--;
}


fdg_loop :: fdg_loop(fdg_settings *settings)
{
  (*this).settings = settings;
  state         = (double*)malloc(2*(*settings).n*sizeof(double));
  energy        = (double*)malloc((*settings).n*sizeof(double));
  delta_energy  = (double*)malloc((*settings).n*sizeof(double));
  distances     = (double*)malloc((((*settings).n*((*settings).n-1))/2)*sizeof(double));
  new_distances = (double*)malloc((*settings).n*sizeof(double));
  num = (*settings).n;

  for (int i = 0; i < (*settings).n; i++){
    state[2*i] = rand_double();
    state[2*i + 1] = rand_double();
    energy[i] = 0;
  }
}

fdg_loop :: ~fdg_loop()
{
  free(state);
  free(energy);
  free(delta_energy);
  free(distances);
  free(new_distances);
}


int fdg_loop :: get_index(int i,int j){
    //indexing into the distance array
    if (i <= j){
      return (((2*(*settings).n-i+1)*i)/2 + (j-i));
    }
    return (((2*(*settings).n-j+1)*j)/2 + (i-j));
  }


double fdg_loop :: update(double beta)
{
    index = 1 + rand()%((*settings).n-1);
    old_x = state[2*index];
    old_y = state[2*index+1];
    new_x = old_x + 2*(rand_double()-0.5)*(*settings).epsilon;
    new_y = old_y + 2*(rand_double()-0.5)*(*settings).epsilon;
    torus_mod(&new_x,&new_y);

    old_energy = 0; new_energy = 0;
    total_delta_energy = 0;
    for (int i = 0; i < (*settings).n; i++){
      new_distances[i] = torus_distance(new_x,new_y,state[2*i],state[2*i+1]);
      delta_energy[i] = (positive_modulo(i-index,(*settings).n) == 1 ? (1-(*settings).A) : 1)*(new_distances[i]-distances[get_index(i,index)]);
      if (i != index){
        total_delta_energy += delta_energy[i];
      }
    }
    if (rand_double() < exp(-beta*total_delta_energy)){
      state[2*index] = new_x;
      state[2*index+1] = new_y;
      for (int i = 0; i < (*settings).n; i++){
        distances[get_index(i,index)] = new_distances[i];
      }
      return total_delta_energy;
    }
    return 0;
}







/*
vector<double[2]> read_boundary(std::string boundary_path){
  //Stored as a list of points in counterclockwise order, e.g. a square is
  //0 0\n1 0\n1 1\n0 1\n

  std::string read_line;
  ifstream boundary_file (boundary_path);

  if (boundary_file.is_open()){
    while (getline(myfile,read_line)){
      strtok(read_line," ");
    }
    boundary_file.close();
  }
}
*/
