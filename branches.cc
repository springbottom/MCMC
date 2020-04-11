
#include <iostream>
#include <string>
#include <fstream>
#include <vector>

int main(int argc; char *argv[]){

  //I'll work with toroidal boundary conditions for now.
  //std::string boundary_path; //The file name for where the boundary is stored.

  MCMC_settings settings;

  if (argc == 7){
    settings.n = argv[1];
    settings.pre_iters = argv[2];
    settings.iters = argv[3];
    settings.temp_min = argv[4];
    settings.temp_max = argv[5];
    settings.temp_N = argv[6];
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
