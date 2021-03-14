#include<stdio.h>
#include <gmsh.h>

int main()
{
  gmsh::initialize();
  gmsh::merge("xpeHb.stl");
  //we are lucky: the model gets geometrized automatically with _1_ surface generated
  gmsh::model::geo::addSurfaceLoop({1},1);
  gmsh::model::geo::addVolume({1},1);
  gmsh::model::geo::synchronize();
  gmsh::model::mesh::generate(3);
  gmsh::fltk::run();
  gmsh::finalize();
  return 0;
}

