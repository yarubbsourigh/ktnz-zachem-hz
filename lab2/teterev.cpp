#include <stdio.h>
#include <cmath>
#include <vector>

#include <vtk-7.1/vtkDoubleArray.h>
#include <vtk-7.1/vtkPoints.h>
#include <vtk-7.1/vtkPointData.h>
#include <vtk-7.1/vtkTetra.h>
#include <vtk-7.1/vtkXMLUnstructuredGridWriter.h>
#include <vtk-7.1/vtkUnstructuredGrid.h>
#include <vtk-7.1/vtkSmartPointer.h>

#include <gmsh.h>

using namespace std;


std::vector<double> vfield(double x, double y, double z, double t) {
  x=x-5;
  y=y+10;
  z=z-4;
  return {-0.1*y+0.1*(x*cos(0.1*t)+y*sin(0.1*t))*cos(0.1*t)*cos(0.2*t),0.1*x+0.1*(x*cos(0.1*t)+y*sin(0.1*t))*sin(0.1*t)*cos(0.2*t),0.1*cos(0.7*x+0.7*y+t)*sin(0.7*x-0.7*y+t)*z};
}

double sfield(double x, double y, double z, double t) {
  x=x-5;
  x=x*6;
  y=y+10;
  y=y*6;
  z=z-4;
  z=z*6;
  return 5*exp(-(pow(sin((0.1*y-0.1*t)/3),2)+pow(sin((0.1*x+0.1*y-0.17*t)/3),2)+pow(sin((0.1*x+0.1*y+0.1*z-0.5*t)/3),2)));
}

std::vector<double> next(double x, double y, double z, double t, double dlt) { //Euler is nothing; Runge-Kutta is everything
  std::vector<double> dr1=vfield(x,y,z,t);
  std::vector<double> dr2=vfield(x+dlt*dr1[0]/2,y+dlt*dr1[1]/2,z+dlt*dr1[2]/2,t+dlt/2);
  std::vector<double> dr3=vfield(x+dlt*dr2[0]/2,y+dlt*dr2[1]/2,z+dlt*dr2[2]/2,t+dlt/2);
  std::vector<double> dr4=vfield(x+dlt*dr3[0],y+dlt*dr3[1],z+dlt*dr3[2],t+dlt);
  return {x+dlt*(dr1[0]+2*dr2[0]+2*dr3[0]+dr4[0])/6,y+dlt*(dr1[1]+2*dr2[1]+2*dr3[1]+dr4[1])/6,z+dlt*(dr1[2]+2*dr2[2]+2*dr3[2]+dr4[2])/6};
}

int main()
{

    double tau = 1;


    gmsh::initialize();
    gmsh::model::add("xpeHb");
    gmsh::merge("xpeHb.stl");
    gmsh::model::geo::addSurfaceLoop({1},1);
    gmsh::model::geo::addVolume({1},1);
    gmsh::model::geo::synchronize();
    gmsh::model::mesh::generate(3);

    vtkUnstructuredGrid *canvas = vtkUnstructuredGrid::New();
    vtkPoints *pts = vtkPoints::New();
    vtkDoubleArray *fnct = vtkDoubleArray::New();
    fnct->SetName("Pa3BoDb1");
    vtkDoubleArray *vel = vtkDoubleArray::New();
    vel->SetName("velocity");
    vel->SetNumberOfComponents(3);

    std::vector<double> ptsxyz;
    std::vector<std::size_t> ptsidx;
    std::vector<double> _;
    gmsh::model::mesh::getNodes(ptsidx, ptsxyz, _, -1, -1, true, false);
    int N = ptsidx.size();

    std::vector<std::size_t> tetidx;
    std::vector<std::size_t> tetvrt;
    gmsh::model::mesh::getElementsByType(4, tetidx, tetvrt);
    int M=tetidx.size();

    gmsh::finalize();

    for(int i = 0; i < N; i++) {
      pts->InsertNextPoint(ptsxyz[3*i], ptsxyz[3*i+1], ptsxyz[3*i+2]);
      double *v = &(vfield(ptsxyz[3*i], ptsxyz[3*i+1], ptsxyz[3*i+2],0)[0]);
      vel->InsertNextTuple(&(vfield(ptsxyz[3*i], ptsxyz[3*i+1], ptsxyz[3*i+2],0)[0]));
      fnct->InsertNextValue(sfield(ptsxyz[3*i], ptsxyz[3*i+1], ptsxyz[3*i+2],0));
    }
    canvas->SetPoints(pts);
    canvas->GetPointData()->AddArray(vel);
    canvas->GetPointData()->AddArray(fnct);

    for(int i=0;i<M; i++) {
      vtkIdList *a = vtkIdList::New();
      for (int j=0; j<4; j++) {
	a->InsertId(j,tetvrt[4*i+j]-1);
      }
      canvas->InsertNextCell(10,a);
    }

    vtkXMLUnstructuredGridWriter *push_backin = vtkXMLUnstructuredGridWriter::New();
    push_backin->SetInputData(canvas);


    for (int i=0; i<180; i=i+1) {
      printf("Generating shot #%d...",i);
      for (int j=0; j<N; j=j+1) {
	double *r = pts->GetPoint(j);
	vector<double> r_ = next(r[0],r[1],r[2],i*tau,tau);
	pts->SetPoint(j,r_[0],r_[1],r_[2]);
	fnct->SetValue(j,sfield(r_[0],r_[1],r_[2],i*tau));
	vel->SetTuple(j,&(vfield(r_[0],r_[1],r_[2],i*tau)[0]));
      }
      printf("Done\n");
      push_backin->SetFileName(("xpeHb-"+std::to_string(i)+".vtu").c_str());
      push_backin->Write();
			    
    }


    return 0;
}
