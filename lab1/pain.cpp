#include <gmsh.h>

int main () {
  gmsh::initialize();
  int centralcenter=1;
  int frontcenter=2;
  int backcenter=3;
  gmsh::model::geo::addPoint(0,0,0,centralcenter);
  gmsh::model::geo::addPoint(10,0,0,frontcenter);
  gmsh::model::geo::addPoint(-10,0,0,backcenter);
  gmsh::model::geo::addPoint(15,0,0,0.5,4);
  gmsh::model::geo::addPoint(5,0,0,0.5,5);
  gmsh::model::geo::addPoint(-15,0,0,0.5,6);
  gmsh::model::geo::addPoint(-5,0,0,0.5,7);
  gmsh::model::geo::addCircleArc(4,centralcenter,6,1);
  gmsh::model::geo::addCircleArc(6,centralcenter,4,2);
  gmsh::model::geo::addCircleArc(5,centralcenter,7,3);
  gmsh::model::geo::addCircleArc(7,centralcenter,5,4);
  gmsh::model::geo::addCircleArc(4,frontcenter,5,5,0,1,0);
  gmsh::model::geo::addCircleArc(5,frontcenter,4,6,0,1,0);
  gmsh::model::geo::addCircleArc(6,backcenter,7,7,0,1,0);
  gmsh::model::geo::addCircleArc(7,backcenter,6,8,0,1,0);
  gmsh::model::geo::addCurveLoop({1,7,-3,6},1);
  gmsh::model::geo::addSurfaceFilling({1},1);
  gmsh::model::geo::addCurveLoop({-2,7,4,6},2);
  gmsh::model::geo::addSurfaceFilling({2},2);
  gmsh::model::geo::addCurveLoop({1,-8,-3,-5},3);
  gmsh::model::geo::addSurfaceFilling({3},3);
  gmsh::model::geo::addCurveLoop({-2,-8,4,-5},4);
  gmsh::model::geo::addSurfaceFilling({4},4);
  gmsh::model::geo::addSurfaceLoop({1,-2,-3,4},1);

  gmsh::model::geo::addPoint(12,0,0,0.5,8);
  gmsh::model::geo::addPoint(8,0,0,0.5,9);
  gmsh::model::geo::addPoint(-12,0,0,0.5,10);
  gmsh::model::geo::addPoint(-8,0,0,0.5,11);
  gmsh::model::geo::addCircleArc(8,centralcenter,10,9);
  gmsh::model::geo::addCircleArc(10,centralcenter,8,10);
  gmsh::model::geo::addCircleArc(9,centralcenter,11,11);
  gmsh::model::geo::addCircleArc(11,centralcenter,9,12);
  gmsh::model::geo::addCircleArc(8,frontcenter,9,13,0,1,0);
  gmsh::model::geo::addCircleArc(9,frontcenter,8,14,0,1,0);
  gmsh::model::geo::addCircleArc(10,backcenter,11,15,0,1,0);
  gmsh::model::geo::addCircleArc(11,backcenter,10,16,0,1,0);
  gmsh::model::geo::addCurveLoop({9,15,-11,14},5);
  gmsh::model::geo::addSurfaceFilling({5},5);
  gmsh::model::geo::addCurveLoop({-10,15,12,14},6);
  gmsh::model::geo::addSurfaceFilling({6},6);
  gmsh::model::geo::addCurveLoop({9,-16,-11,-13},7);
  gmsh::model::geo::addSurfaceFilling({7},7);
  gmsh::model::geo::addCurveLoop({-10,-16,12,-13},8);
  gmsh::model::geo::addSurfaceFilling({8},8);
  gmsh::model::geo::addSurfaceLoop({5,-6,-7,8},2);
  gmsh::model::geo::addVolume({1,-2},1);
  gmsh::model::geo::synchronize();
  gmsh::model::mesh::generate(3);
  gmsh::fltk::run();
  gmsh::finalize();
  return 0;

}
