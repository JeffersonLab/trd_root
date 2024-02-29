#include <stdio.h>
#include <string>
#include <vector>

void printVector(std::string name, std::vector<int>& vec){
  printf("1D Vector: %s\n",name.c_str());
  int cnt=0;
  for(int x = 0; x < vec.size(); x++){
    printf("i=%d %d \n ",cnt,vec[x]); cnt++;
  }
  printf("\n");
}

void printVector(std::string name, std::vector<float>& vec){
  printf("1D Vector: %s\n",name.c_str());
  int cnt=0;
  for(int x = 0; x < vec.size(); x++){
    printf("k=%d %3.0f, ",cnt,vec[x]); cnt++;
  }
  printf("\n");
}

void printVector(std::string name, std::vector<std::vector<float>>& vec){
  printf("2D Vector: %s\n",name.c_str());
  int cnt=0;
  for(int x = 0; x < vec.size(); x++){
    printf(" n=%d  ",cnt); 
    for(int y = 0; y < vec[x].size(); y++){
      printf(" %6.3f ",vec[x][y]); 
    }
    printf("\n"); cnt++;
  }
  printf("\n");
}

void printVector(std::string name, std::vector<std::vector<std::vector<float>>>& vec){
  printf("3D Vector: %s\n",name.c_str());
  for(int x = 0; x < vec.size(); x++){
    for(int y = 0; y < vec[x].size(); y++){
      for(int z = 0; z < 3; z++){
        printf(" %6.3f ",vec[x][y][z]);
      }
      printf("\n");
    }
    printf("|\n");
  }
  printf("\n");
}


void printToPython(std::vector<std::vector<float>> v){
  std::string printFormat = "%5.3f";
  printf("---------------\n");
  int max = 0;
  for(int i = 0; i < v.size(); i++){
    max = std::max(int(v[i].size()),max);
  }

  for(int i = 0; i < v.size(); i++){
    for(int j = 0; j < v[i].size(); j++){
      if(j < max-1){
        printf((printFormat+",").c_str(),v[i][j]);
      } else {
        printf(printFormat.c_str(),v[i][j]);
      }
    }
    for(int j = v[i].size(); j < max; j++){
      if(j < max-1){
        printf((printFormat+",").c_str(),0);
      } else {
        printf(printFormat.c_str(),0);
      }
    }
    printf("\n");
  }
  printf("---------------\n");
}

void printToPython(std::vector<std::vector<int>> v){
  std::string printFormat = "%3d";
  printf("---------------\n");
  int max = 0;
  for(int i = 0; i < v.size(); i++){
    max = std::max(int(v[i].size()),max);
  }

  for(int i = 0; i < v.size(); i++){
    for(int j = 0; j < v[i].size(); j++){
      if(j < max-1){
        printf((printFormat+",").c_str(),v[i][j]);
      } else {
        printf(printFormat.c_str(),v[i][j]);
      }
    }
    for(int j = v[i].size(); j < max; j++){
      if(j < max-1){
        printf((printFormat+",").c_str(),-1);
      } else {
        printf(printFormat.c_str(),-1);
      }
    }
    printf("\n");
  }
  printf("---------------\n");
}
