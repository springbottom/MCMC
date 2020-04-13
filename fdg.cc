#include "MCMC.h"
#include "fdg.h"

int main(int argc, char *argv[]){

  MCMC_settings settings;
  fdg_settings settings2;

  if (argc == 7){
    settings2.n = atoi(argv[1]);
    settings2.epsilon = strtod(argv[2],NULL);
    settings2.A = strtod(argv[3],NULL);

    settings.pre_iters = atoi(argv[4]);
    settings.iters = atoi(argv[5]);
    settings.temp_min = strtod(argv[6],NULL);
    settings.temp_max = strtod(argv[7],NULL);
    settings.temp_N = atoi(argv[8]);
  }

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

class fdg_loop{
  double* state;
  double* energy; //represents all the bonds connected to that vertex...
  double* delta_energy;
  double* distances;
  fdg_settings* settings;

  fdg_loop(){
    state  = (double*)malloc(2*(*settings).n*sizeof(double));
    energy = (double*)malloc((*settings).n*sizeof(double));
    distances = (double*)malloc((*settings).n*(*settings).n*sizeof(double));
    new_distances = (double*)malloc((*settings).n*sizeof(double));
  }

  int index;
  double new_x, new_y;
  double old_x, old_y;
  double* new_distances;
  double old_energy, new_energy;
  double total_delta_energy;

  int get_index(int i,int j){
    //indexing into the distance array
    if (i <= j){
      return (2*(*settings).n-i+1)*(i/(double)2) + j;
    }
    return (2*(*settings).n-j+1)*(j/(double)2) + i;
  }


  double update(double beta){
    index = 1 + rand()%((*settings).n-1);
    old_x = state[2*index];
    old_y = state[2*index+1];
    new_x = old_x + 2*(rand_double()-0.5)*(*settings).epsilon;
    new_y = old_y + 2*(rand_double()-0.5)*(*settings).epsilon;
    new_x = (new_x < 0 ? new_x + 1 : new_x);
    new_x = (new_x > 1 ? new_x - 1 : new_x);
    new_y = (new_y < 0 ? new_y + 1 : new_y);
    new_y = (new_y > 1 ? new_y - 1 : new_y);

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

  int init(){
    for (int i = 0; i < (*settings).n; i++){
      state[2*i] = rand_double();
      state[2*i + 1] = rand_double();
      energy[i] = 0;
    }
    return 0;
  }
};






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
