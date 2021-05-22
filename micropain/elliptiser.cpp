//compile with -lsfml-graphics -lsfml-window -lsfml-system
//Enjoy :)
#include <stdio.h>
#include <math.h>
#include<string.h>
#include<vector>
using namespace std;
#include<SFML/Graphics.hpp>

double pi = 3.1415926535897932384626433832795;

typedef struct point {
  double x;
  double y;
} point;

typedef struct arrow {
  double x;
  double y;
} arrow;

typedef struct form {
    double xx;
    double xy;
    double yy;
  } form;

typedef struct ellipse {
  form shape;
  point center;
  double threshold;
} ellipse;

arrow fromto(point A, point B) {
  arrow a;
  a.x=B.x-A.x;
  a.y=B.y-A.y;
  return a;
}

double area(ellipse E) {
  return pi*E.threshold/sqrt(E.shape.xx*E.shape.yy-E.shape.xy*E.shape.xy);
}

double evl(form F, arrow p, arrow q) {
  return F.xx*p.x*q.x+F.xy*p.x*q.y+F.xy*p.y*q.x+F.yy*p.y*q.y;
}

double cross(arrow A, arrow B) {
  return A.x*B.y-A.y*B.x;
}

bool tudee(point A, point B) {
  if (A.x != B.x)
    return A.x>B.x;
  else return A.y>B.y;
}

sf::VertexArray *drellipse(ellipse E,int N) {
  sf::VertexArray *oel = new sf::VertexArray(sf::LineStrip,N+1);
  for (int i=0; i<N; i=i+1) {
    double phi = 2*pi*i/N;
    arrow tst; tst.x = cos(phi); tst.y=sin(phi);
    double NF= sqrt(E.threshold/evl(E.shape,tst,tst));
    tst.x=NF*tst.x; tst.y=NF*tst.y;
    (*oel)[i].position = sf::Vector2f(E.center.x+tst.x,-(E.center.y+tst.y));
    (*oel)[i].color=sf::Color::Blue;
  }
  (*oel)[N].position=(*oel)[0].position;
  (*oel)[N].color=sf::Color::Blue;
  return oel;
}

ellipse threep(point A, point B, point C) {
  point centr;
  centr.x = (A.x+B.x+C.x)/3;
  centr.y = (A.y+B.y+C.y)/3;
  arrow p = fromto(centr,A);
  arrow q = fromto(centr,B);
  arrow r = fromto(centr,C);
  form F;
  F.xx = p.y*p.y+q.y*q.y+r.y*r.y;
  F.xy = -p.x*p.y-q.x*q.y-r.x*r.y;
  F.yy = p.x*p.x+q.x*q.x+r.x*r.x;
  double thr = evl(F,p,p);
  ellipse E;
  E.shape=F;
  E.center=centr;
  E.threshold=thr;
  return E;
}

ellipse fourp(point A, point B, point C, point D) {
  double a = cross(fromto(A,C),fromto(A,B));
  double b = -cross(fromto(A,C),fromto(A,D));
  double c = -cross(fromto(B,D),fromto(B,A));
  double d = cross(fromto(B,D),fromto(B,C));
  point dia;
  dia.x = (B.x*b+D.x*a)/(a+b);
  dia.y = (B.y*b+D.y*a)/(a+b);
  
  a=a/sqrt(a*b);
  b=1/a;
  c=c/sqrt(c*d);
  d=1/c;
  double r = pi*0.99999;
  double l = 0.00001;
  while (r-l>1e-8) {
    double x = (r+l)/2;
    if ((2*a*b*cos(x)-(a-b)*(d-c)/2)*sin(x)*sin(x)>3*cos(x)*(a*b*sin(x)*sin(x)+((a-b)*(a-b)/4)+((d-c)*(d-c)/4)+(a-b)*(d-c)*cos(x)/2))
      r=x;
    else
      l=x;
  }
  double x=(l+r)/2;
  arrow f; f.x=0; f.y=a;
  arrow g; g.x=d*sin(x); g.y=-d*cos(x);
  arrow F = fromto(dia,B);
  arrow G = fromto(dia,C);
  arrow S;
  S.x = ((a-b)/2)*cos(x)/sin(x)+((d-c)/2)*1/sin(x);
  S.y = ((a-b)/2);
  double Sf = (g.y*S.x-g.x*S.y)/(f.x*g.y-f.y*g.x);
  double Sg = (-f.y*S.x+f.x*S.y)/(f.x*g.y-f.y*g.x);
  point O;
  O.x = dia.x+Sf*F.x+Sg*G.x;
  O.y = dia.y+Sf*F.y+Sg*G.y;
  
  form M;
  M.xx = G.y*G.y*a*a+F.y*F.y*d*d+2*F.y*G.y*a*d*cos(x);
  M.yy = G.x*G.x*a*a+F.x*F.x*d*d+2*F.x*G.x*a*d*cos(x);
  M.xy = -G.y*G.x*a*a-F.y*F.x*d*d-(F.x*G.y+F.y*G.x)*a*d*cos(x);

  ellipse E;
  E.shape = M;
  E.center = O;
  E.threshold = evl(M,fromto(O,A),fromto(O,A));

  return E;
}

ellipse fivep(point A, point B, point C, point D, point E) {
  double x_1,x_2,y_1,y_2,z_1,z_2,u_1,u_2,v_1,v_2;
  x_1=A.x; x_2=A.y;
  y_1=B.x; y_2=B.y;
  z_1=C.x; z_2=C.y;
  u_1=D.x; u_2=D.y;
  v_1=E.x; v_2=E.y;
  double detxx, detxy, detyy, detx, dety,det1;
  
  detxx=-v_2*x_2*y_2*(v_2*(x_1-y_1)+x_2*(-v_1+y_1)+(v_1-x_1)*y_2)*(u_1-z_1)+(v_2*(x_2*x_2*y_1-x_1*y_2*y_2)*(v_1-z_1)+v_1*x_2*y_2*(y_2*(x_1-z_1)+x_2*(-y_1+z_1))+v_2*v_2*(x_2*y_1*z_1-x_1*(y_1*(x_2-y_2)+y_2*z_1))+u_1*(-v_2*(x_2*x_2-y_2*y_2)*(v_1-z_1)+v_2*v_2*(x_1*x_2-y_1*y_2+(-x_2+y_2)*z_1)+x_2*y_2*(x_2*(y_1-z_1)+y_2*(-x_1+z_1))))*z_2+(v_2*(v_1-x_1)*x_2*(u_1-y_1)+(-v_2*x_1*y_1+u_1*(-v_1*v_2+x_2*(x_1-y_1)+v_2*y_1)+v_1*(x_1*(v_2-x_2)+x_2*y_1))*y_2)*z_2*z_2+u_2*u_2*(x_2*(-x_1+y_1)*y_2*z_1+(-x_2*y_1*z_1+x_1*(y_1*(x_2-y_2)+y_2*z_1))*z_2+v_1*(x_2*(x_1-y_1)*y_2+(-x_1*x_2+y_2*(y_1-z_1)+x_2*z_1)*z_2+v_2*(x_2*(y_1-z_1)+y_2*(-x_1+z_1)+(x_1-y_1)*z_2))+v_2*(y_1*z_1*(-y_2+z_2)+x_1*(y_1*y_2+x_2*(-y_1+z_1)-z_1*z_2)))+u_2*(x_2*y_2*(x_2*y_1-x_1*y_2)*(v_1-z_1)+(x_2*x_2*(-v_1+y_1)+(v_1-x_1)*y_2*y_2)*z_1*z_2+(x_1*x_2*(v_1-y_1)+(-v_1+x_1)*y_1*y_2)*z_2*z_2+v_1*v_2*(y_2*y_2*(x_1-z_1)+x_2*x_2*(-y_1+z_1)+(-x_1+y_1)*z_2*z_2)+v_2*v_2*(y_1*z_1*(y_2-z_2)+x_1*(-y_1*y_2+x_2*(y_1-z_1)+z_1*z_2))+u_1*(x_2*y_2*(-x_2+y_2)*(v_1-z_1)+(x_2*x_2*(v_1-y_1)+(-v_1+x_1)*y_2*y_2)*z_2+(x_2*(-v_1+y_1)+(v_1-x_1)*y_2)*z_2*z_2+v_2*v_2*(y_2*(x_1-z_1)+x_2*(-y_1+z_1)+(-x_1+y_1)*z_2)+v_2*(x_2*x_2*(y_1-z_1)+y_2*y_2*(-x_1+z_1)+(x_1-y_1)*z_2*z_2)));
  
  detxy = z_1*(v_1*x_2*(x_2-y_2)*y_2*(v_1-z_1)+v_2*v_2*(x_2*y_1*(y_1-z_1)+x_1*y_2*(-x_1+z_1))+v_2*(x_1*y_2*y_2*(x_1-z_1)+x_2*x_2*y_1*(-y_1+z_1)))+(v_2*v_2*x_1*(x_1-y_1)*y_1+v_1*(x_2*x_2*y_1*(-v_1+y_1)+(v_1-x_1)*x_1*y_2*y_2))*z_2+(v_2*x_1*y_1*(-x_1+y_1)+v_1*(x_2*(v_1-y_1)*y_1+x_1*(-v_1+x_1)*y_2))*z_2*z_2+u_2*u_2*(v_2*(x_1-y_1)*(x_1-z_1)*(y_1-z_1)+(v_1-z_1)*(x_2*y_1*(y_1-z_1)+x_1*y_2*(-x_1+z_1)+v_1*(y_2*(x_1-z_1)+x_2*(-y_1+z_1)))+(v_1-x_1)*(v_1-y_1)*(-x_1+y_1)*z_2)+u_1*u_1*(x_2*(x_2-y_2)*y_2*(v_1-z_1)+(x_2*x_2*(-v_1+y_1)+(v_1-x_1)*y_2*y_2)*z_2+(x_2*(v_1-y_1)+(-v_1+x_1)*y_2)*z_2*z_2+v_2*v_2*(x_2*(y_1-z_1)+y_2*(-x_1+z_1)+(x_1-y_1)*z_2)+v_2*(y_2*y_2*(x_1-z_1)+x_2*x_2*(-y_1+z_1)+(-x_1+y_1)*z_2*z_2))+u_1*(x_2*(x_2-y_2)*y_2*z_1*z_1-v_1*v_1*(x_2-y_2)*(x_2-z_2)*(y_2-z_2)+(-x_2*x_2*y_1*y_1+x_1*x_1*y_2*y_2)*z_2+(x_2*y_1*y_1-x_1*x_1*y_2)*z_2*z_2+v_2*v_2*(y_2*(x_1*x_1-z_1*z_1)+x_2*(-y_1*y_1+z_1*z_1)+(-x_1*x_1+y_1*y_1)*z_2)+v_2*(x_2*x_2*(y_1*y_1-z_1*z_1)+y_2*y_2*(-x_1*x_1+z_1*z_1)+(x_1*x_1-y_1*y_1)*z_2*z_2))+u_2*(-v_2*v_2*(x_1-y_1)*(x_1-z_1)*(y_1-z_1)+z_1*(x_2*x_2*y_1*(y_1-z_1)+x_1*y_2*y_2*(-x_1+z_1))+x_1*(x_1-y_1)*y_1*z_2*z_2+v_1*v_1*(x_2*x_2*(y_1-z_1)+y_2*y_2*(-x_1+z_1)+(x_1-y_1)*z_2*z_2)+v_1*(y_2*y_2*(x_1*x_1-z_1*z_1)+x_2*x_2*(-y_1*y_1+z_1*z_1)+(-x_1*x_1+y_1*y_1)*z_2*z_2));

    detyy =  z_1*(v_1*v_1*x_2*(-x_1+y_1)*y_2+v_2*x_1*y_1*(x_2*(y_1-z_1)+y_2*(-x_1+z_1))+v_1*(x_2*(x_1-y_1)*y_2*z_1+v_2*(x_1*y_2*(x_1-z_1)+x_2*y_1*(-y_1+z_1))))+(v_1*x_1*y_1*(x_2*(v_1-y_1)+v_2*(-x_1+y_1)+(-v_1+x_1)*y_2)+(v_2*x_1*(x_1-y_1)*y_1+v_1*(x_2*y_1*(-v_1+y_1)+(v_1-x_1)*x_1*y_2))*z_1)*z_2+u_1*(-v_2*x_1*x_2*y_1*y_1+v_2*x_1*x_1*y_1*y_2+v_2*x_1*x_2*z_1*z_1-x_1*x_2*y_2*z_1*z_1-v_2*y_1*y_2*z_1*z_1+x_2*y_1*y_2*z_1*z_1+(x_1*x_2*y_1*y_1+(v_2-x_2)*y_1*y_1*z_1-x_1*x_1*(y_1*y_2+(v_2-y_2)*z_1))*z_2+u_2*(-v_2*(x_1-y_1)*(x_1-z_1)*(y_1-z_1)+(v_1-z_1)*(x_1*y_2*(x_1-z_1)+x_2*y_1*(-y_1+z_1)+v_1*(x_2*(y_1-z_1)+y_2*(-x_1+z_1)))+(v_1-x_1)*(v_1-y_1)*(x_1-y_1)*z_2)+v_1*v_1*(x_2*(x_1-y_1)*y_2+(-x_1*x_2+y_1*y_2+(x_2-y_2)*z_1)*z_2)+v_1*v_2*(y_2*z_1*z_1+x_2*(y_1*y_1-z_1*z_1)-y_1*y_1*z_2+x_1*x_1*(-y_2+z_2)))+u_2*(x_1*y_1*z_1*(y_2*(x_1-z_1)+x_2*(-y_1+z_1)+(-x_1+y_1)*z_2)+v_1*v_1*(y_1*z_1*(-y_2+z_2)+x_1*(y_1*y_2+x_2*(-y_1+z_1)-z_1*z_2))+v_1*(v_2*(x_1-y_1)*(x_1-z_1)*(y_1-z_1)+x_1*x_2*(y_1*y_1-z_1*z_1)+y_1*z_1*(y_2*z_1-y_1*z_2)+x_1*x_1*(-y_1*y_2+z_1*z_2)))+u_1*u_1*(x_2*(x_1-y_1)*y_2*z_1+(x_2*y_1*z_1-x_1*(y_1*(x_2-y_2)+y_2*z_1))*z_2+v_1*(x_2*(-x_1+y_1)*y_2+(x_1*x_2-y_1*y_2+(-x_2+y_2)*z_1)*z_2+v_2*(y_2*(x_1-z_1)+x_2*(-y_1+z_1)+(-x_1+y_1)*z_2))+v_2*(y_1*z_1*(y_2-z_2)+x_1*(-y_1*y_2+x_2*(y_1-z_1)+z_1*z_2)));

    detx =  v_2*x_2*y_2*(v_2*(x_1-y_1)+x_2*(-v_1+y_1)+(v_1-x_1)*y_2)*(u_1*u_1-z_1*z_1)+((v_2*x_1-v_1*x_2)*(v_2*y_1-v_1*y_2)*(x_2*y_1-x_1*y_2)+(v_1*v_1*x_2*y_2*(-x_2+y_2)+v_2*v_2*(-x_2*y_1*y_1+x_1*x_1*y_2)+v_2*(x_2*x_2*y_1*y_1-x_1*x_1*y_2*y_2))*z_1+u_1*u_1*(v_2*(x_2*x_2-y_2*y_2)*(v_1-z_1)+v_2*v_2*(-x_1*x_2+y_2*(y_1-z_1)+x_2*z_1)+x_2*y_2*(y_2*(x_1-z_1)+x_2*(-y_1+z_1))))*z_2+(-v_2*(v_1-x_1)*x_2*(u_1*u_1-y_1*y_1)+(-v_1*v_2*x_1*x_1+v_1*v_1*x_2*(x_1-y_1)+v_2*x_1*x_1*y_1+u_1*u_1*(v_1*v_2-x_1*x_2+(-v_2+x_2)*y_1))*y_2)*z_2*z_2+u_2*u_2*(x_2*(x_1-y_1)*y_2*z_1*z_1+(-x_1*x_2*y_1*y_1+x_1*x_1*y_2*(y_1-z_1)+x_2*y_1*y_1*z_1)*z_2+v_1*v_2*(y_2*(x_1*x_1-z_1*z_1)+x_2*(-y_1*y_1+z_1*z_1)+(-x_1*x_1+y_1*y_1)*z_2)+v_1*v_1*(x_2*(-x_1+y_1)*y_2+(x_1*x_2-y_1*y_2+(-x_2+y_2)*z_1)*z_2)+v_2*(x_1*x_2*(y_1*y_1-z_1*z_1)+y_1*z_1*(y_2*z_1-y_1*z_2)+x_1*x_1*(-y_1*y_2+z_1*z_2)))+u_2*(-x_2*y_2*(x_2*y_1-x_1*y_2)*(v_1*v_1-z_1*z_1)+(x_2*x_2*(v_1*v_1-y_1*y_1)+(-v_1*v_1+x_1*x_1)*y_2*y_2)*z_1*z_2+(x_1*x_2*(-v_1*v_1+y_1*y_1)+(v_1*v_1-x_1*x_1)*y_1*y_2)*z_2*z_2+v_1*v_2*(x_2*x_2*(y_1*y_1-z_1*z_1)+y_2*y_2*(-x_1*x_1+z_1*z_1)+(x_1*x_1-y_1*y_1)*z_2*z_2)+v_2*v_2*(x_1*x_2*(-y_1*y_1+z_1*z_1)+y_1*z_1*(-y_2*z_1+y_1*z_2)+x_1*x_1*(y_1*y_2-z_1*z_2))+u_1*(x_2*y_2*(-x_2+y_2)*z_1*z_1+v_1*v_1*(x_2-y_2)*(x_2-z_2)*(y_2-z_2)+(x_2*x_2*y_1*y_1-x_1*x_1*y_2*y_2)*z_2+(-x_2*y_1*y_1+x_1*x_1*y_2)*z_2*z_2+v_2*v_2*(x_2*(y_1*y_1-z_1*z_1)+y_2*(-x_1*x_1+z_1*z_1)+(x_1*x_1-y_1*y_1)*z_2)+v_2*(y_2*y_2*(x_1*x_1-z_1*z_1)+x_2*x_2*(-y_1*y_1+z_1*z_1)+(-x_1*x_1+y_1*y_1)*z_2*z_2)));


    dety = z_1*(-v_1*x_2*y_2*(x_2*y_1-x_1*y_2)*(v_1-z_1)+v_1*v_2*(x_2*x_2*y_1*(y_1-z_1)+x_1*y_2*y_2*(-x_1+z_1))+v_2*v_2*x_1*y_1*(y_2*(x_1-z_1)+x_2*(-y_1+z_1)))+(v_2*v_2*x_1*y_1*(-x_1+y_1)+v_1*(x_2*x_2*(v_1-y_1)*y_1+x_1*(-v_1+x_1)*y_2*y_2))*z_1*z_2+v_1*x_1*y_1*(v_2*(x_1-y_1)+x_2*(-v_1+y_1)+(v_1-x_1)*y_2)*z_2*z_2+u_1*u_2*(v_2*v_2*(x_1-y_1)*(x_1-z_1)*(y_1-z_1)+z_1*(x_1*y_2*y_2*(x_1-z_1)+x_2*x_2*y_1*(-y_1+z_1))+x_1*y_1*(-x_1+y_1)*z_2*z_2+v_1*v_1*(y_2*y_2*(x_1-z_1)+x_2*x_2*(-y_1+z_1)+(-x_1+y_1)*z_2*z_2)+v_1*(x_2*x_2*(y_1*y_1-z_1*z_1)+y_2*y_2*(-x_1*x_1+z_1*z_1)+(x_1*x_1-y_1*y_1)*z_2*z_2))+u_1*u_1*(-x_2*y_2*(x_2*y_1-x_1*y_2)*(v_1-z_1)+(x_2*x_2*(v_1-y_1)+(-v_1+x_1)*y_2*y_2)*z_1*z_2+(x_1*x_2*(-v_1+y_1)+(v_1-x_1)*y_1*y_2)*z_2*z_2+v_1*v_2*(x_2*x_2*(y_1-z_1)+y_2*y_2*(-x_1+z_1)+(x_1-y_1)*z_2*z_2)+v_2*v_2*(y_1*z_1*(-y_2+z_2)+x_1*(y_1*y_2+x_2*(-y_1+z_1)-z_1*z_2)))+u_1*((x_2*y_1-x_1*y_2)*(x_2*z_1-x_1*z_2)*(-y_2*z_1+y_1*z_2)+v_1*v_2*(y_2*y_2*(x_1*x_1-z_1*z_1)+x_2*x_2*(-y_1*y_1+z_1*z_1)+(-x_1*x_1+y_1*y_1)*z_2*z_2)+v_1*v_1*(x_2*y_2*(x_2*y_1-x_1*y_2)+(-x_2*x_2+y_2*y_2)*z_1*z_2+(x_1*x_2-y_1*y_2)*z_2*z_2)+v_2*v_2*(x_1*x_2*(y_1*y_1-z_1*z_1)+y_1*z_1*(y_2*z_1-y_1*z_2)+x_1*x_1*(-y_1*y_2+z_1*z_2)))+u_2*u_2*(x_1*y_1*z_1*(x_2*(y_1-z_1)+y_2*(-x_1+z_1)+(x_1-y_1)*z_2)+v_1*(-v_2*(x_1-y_1)*(x_1-z_1)*(y_1-z_1)+x_1*x_2*(-y_1*y_1+z_1*z_1)+y_1*z_1*(-y_2*z_1+y_1*z_2)+x_1*x_1*(y_1*y_2-z_1*z_2))+v_1*v_1*(y_1*z_1*(y_2-z_2)+x_1*(-y_1*y_2+x_2*(y_1-z_1)+z_1*z_2)));

    det1 =  u_2*u_2*(z_1*(v_1*v_1*x_2*(x_1-y_1)*y_2+v_2*x_1*y_1*(y_2*(x_1-z_1)+x_2*(-y_1+z_1))+v_1*(x_2*(-x_1+y_1)*y_2*z_1+v_2*(x_2*y_1*(y_1-z_1)+x_1*y_2*(-x_1+z_1))))+(v_1*x_1*y_1*(v_2*(x_1-y_1)+x_2*(-v_1+y_1)+(v_1-x_1)*y_2)+(v_2*x_1*y_1*(-x_1+y_1)+v_1*(x_2*(v_1-y_1)*y_1+x_1*(-v_1+x_1)*y_2))*z_1)*z_2)+u_1*(-v_2*x_2*y_2*(v_2*(x_1-y_1)+x_2*(-v_1+y_1)+(v_1-x_1)*y_2)*(u_1-z_1)*z_1+((v_2*x_1-v_1*x_2)*(v_2*y_1-v_1*y_2)*(-x_2*y_1+x_1*y_2)+(v_1*v_1*x_2*(x_2-y_2)*y_2+v_2*v_2*(x_2*y_1*y_1-x_1*x_1*y_2)+v_2*(-x_2*x_2*y_1*y_1+x_1*x_1*y_2*y_2))*z_1+u_1*(-v_2*(x_2*x_2*y_1-x_1*y_2*y_2)*(v_1-z_1)+v_1*x_2*y_2*(x_2*(y_1-z_1)+y_2*(-x_1+z_1))+v_2*v_2*(-x_2*y_1*z_1+x_1*(y_1*(x_2-y_2)+y_2*z_1))))*z_2+(v_2*(v_1-x_1)*x_2*(u_1-y_1)*y_1+(v_1*x_1*(v_2*(-u_1+x_1)+(u_1-v_1)*x_2)+(v_2*(u_1-x_1)*x_1+v_1*(-u_1+v_1)*x_2)*y_1)*y_2)*z_2*z_2)+u_2*(z_1*(v_1*x_2*y_2*(x_2*y_1-x_1*y_2)*(v_1-z_1)+v_2*v_2*x_1*y_1*(x_2*(y_1-z_1)+y_2*(-x_1+z_1))+v_1*v_2*(x_1*y_2*y_2*(x_1-z_1)+x_2*x_2*y_1*(-y_1+z_1)))+(v_2*v_2*x_1*(x_1-y_1)*y_1+v_1*(x_2*x_2*y_1*(-v_1+y_1)+(v_1-x_1)*x_1*y_2*y_2))*z_1*z_2+v_1*x_1*y_1*(x_2*(v_1-y_1)+v_2*(-x_1+y_1)+(-v_1+x_1)*y_2)*z_2*z_2+u_1*(v_1*(x_2*y_2*(-x_2+y_2)*(v_1-z_1)*z_1+(x_2*x_2*(v_1-y_1)*y_1+x_1*(-v_1+x_1)*y_2*y_2)*z_2+(x_2*y_1*(-v_1+y_1)+(v_1-x_1)*x_1*y_2)*z_2*z_2)+v_2*v_2*(x_2*y_1*z_1*(-y_1+z_1)+x_1*(y_2*(x_1-z_1)*z_1+y_1*(-x_1+y_1)*z_2))+v_2*(x_2*x_2*y_1*(y_1-z_1)*z_1+x_1*(y_2*y_2*z_1*(-x_1+z_1)+(x_1-y_1)*y_1*z_2*z_2))));

    form UF;
    UF.xx=detxx;
    UF.xy=detxy/2;
    UF.yy=detyy;
    double delt = UF.xx*UF.yy - UF.xy*UF.xy;
    arrow Uc;
    Uc.x = - (UF.yy*detx-UF.xy*dety)/(2*delt);
    Uc.y = - (-UF.xy*detx+UF.xx*dety)/(2*delt);
    
    double Ut;
    Ut=-det1+evl(UF,Uc,Uc);

    ellipse U;
    U.shape=UF;
    U.center.x=Uc.x; U.center.y=Uc.y;
    U.threshold=Ut;
    return U;
}


bool circum(ellipse E, vector<point> *pts) {
  for (int i=0; i<pts->size(); i=i+1) {
    arrow del = fromto(E.center,(*pts)[i]);
    if (abs(evl(E.shape,del,del))-abs(E.threshold)>1)
      return false;
  }
  return true;
}

ellipse l4three(vector<point> *pts) {
  int N=pts->size();
  ellipse E; E.threshold=0;
  double Ar;
  for (int i=0; i<N; i=i+1)
    for (int j=i+1; j<N; j=j+1) {
      arrow IJ=fromto((*pts)[i],(*pts)[j]);
      if (j-i>1) {
	int L=i; int R=j; int K=(R+L)/2;
	while (R-L>2) {
	  int nK=(K+1)%N;
	  int pK=(N+K-1)%N;
	  arrow naK=fromto((*pts)[K],(*pts)[nK]);
	  arrow paK=fromto((*pts)[pK],(*pts)[K]);
	  if (cross(IJ,paK)*cross(IJ,naK)<=0)
	    break;
	  else if (cross(IJ,paK)>0)
	    L=K;
	  else
	    R=K;
	  K=(R+L)/2;
	}
	if (E.threshold==0 || abs(cross(IJ,fromto((*pts)[i],(*pts)[K])))/2>Ar) {
	  E=threep((*pts)[i],(*pts)[j],(*pts)[K]);
	  Ar = abs(cross(IJ,fromto((*pts)[i],(*pts)[K])))/2;
	}
      }
      if (N-j+i>1) {
	int L=j; int R=N+i; int K=(R+L)/2;
	while (R-L>2) {
	  int k=(K)%N;
	  int nk=(k+1)%N;
	  int pk=(N+k-1)%N;
	  arrow naK=fromto((*pts)[k],(*pts)[nk]);
	  arrow paK=fromto((*pts)[pk],(*pts)[k]);
	  if (cross(IJ,paK)*cross(IJ,naK)<=0)
	    break;
	  else if (cross(IJ,paK)<0)
	    L=K;
	  else
	    R=K;
	  K=(R+L)/2;
	}
	if (E.threshold==0 || abs(cross(IJ,fromto((*pts)[K%N],(*pts)[i])))/2>Ar) {
	  E=threep((*pts)[i],(*pts)[j],(*pts)[K%N]);
	  Ar = abs(cross(IJ,fromto((*pts)[K%N],(*pts)[i])))/2;
	}
      }
    }
  return E;
}

ellipse l4four(vector<point> *pts) {
  int N=pts->size();
  ellipse E; E.threshold=0;
  E.shape.xx=E.shape.xy=E.shape.yy=0;
  E.center.x=E.center.y=0;
  for (int i=0; i<N; i=i+1)
    for (int j=i+1; j<N; j=j+1)
      for (int k=j+1; k<N; k=k+1)
	for (int l=k+1; l<N; l=l+1) {
	  ellipse U=fourp((*pts)[i],(*pts)[j],(*pts)[k],(*pts)[l]);
	  if (circum(U,pts)) {
	    return U;
	  }
	}
  return E;
}

ellipse l4five(vector<point> *pts) {
  int N=pts->size();
  ellipse E; E.threshold=0;
  E.shape.xx=E.shape.xy=E.shape.yy=0;
  E.center.x=E.center.y=0;
  for (int i=0; i<N; i=i+1)
    for (int j=i+1; j<N; j=j+1)
      for (int k=j+1; k<N; k=k+1)
	for (int l=k+1; l<N; l=l+1)
	  for (int m=l+1; m<N; m=m+1) {
	    ellipse U=fivep((*pts)[i],(*pts)[j],(*pts)[k],(*pts)[l],(*pts)[m]);
	    if ((U.shape.xx*U.shape.yy-U.shape.xy*U.shape.xy)>0)
	      if (circum(U,pts))
		return U;
	  }
  return E;
}

ellipse encloser(vector<point> *pts) {
  ellipse E;
  E=l4three(pts);
  if (circum(E,pts))
    return E;
  E=l4four(pts);
  if (E.threshold>0)
    return E;
  E=l4five(pts);
  return E;
}

void ptsort(vector<point> *pts) {
  for (int i=1; i<(*pts).size(); i=i+1) {
    int nw=i;
    while (nw>0 && tudee((*pts)[nw],(*pts)[(nw+1)/2-1])) {
      point t = (*pts)[(nw+1)/2-1];
      (*pts)[(nw+1)/2-1]=(*pts)[nw];
       (*pts)[nw]=t;
      nw=(nw+1)/2-1;
    }
  }
  for (int i=(*pts).size()-1; i>0; i=i-1) {
    point t =(*pts)[i]; (*pts)[i]=(*pts)[0]; (*pts)[0]=t;
    int nw=0;
    while (1) {
      if (2*nw+1>=i)
        break;
      else if (2*nw+2>=i) {
	if (tudee((*pts)[nw],(*pts)[2*nw+1]))
	   break;
	else {
	  point t=(*pts)[2*nw+1]; (*pts)[2*nw+1]=(*pts)[nw]; (*pts)[nw]=t;
	  break;
	}
      }
      else {
	if (tudee((*pts)[nw],(*pts)[2*nw+1]) && tudee((*pts)[nw],(*pts)[2*nw+2]))
	    break;
	else if (tudee((*pts)[2*nw+1],(*pts)[2*nw+2])) {
	  point t=(*pts)[2*nw+1]; (*pts)[2*nw+1]=(*pts)[nw]; (*pts)[nw]=t;
	  nw=2*nw+1;
	}
	else {
	  point t=(*pts)[2*nw+2]; (*pts)[2*nw+2]=(*pts)[nw]; (*pts)[nw]=t;
	  nw=2*nw+2;
	}
      }
    }
  }
}

vector<point> *hull(vector<point> *pts) {
  vector<point> *H = new vector<point>;
   ptsort(pts);
  (*H).push_back((*pts)[0]);
  (*H).push_back((*pts)[1]);
  for (int i=2; i<(*pts).size(); i=i+1) {
    while ((*H).size()>=2 && cross(fromto((*H)[(*H).size()-2],(*H)[(*H).size()-1]),
				fromto((*H)[(*H).size()-1],(*pts)[i]))>=0)
      (*H).pop_back();
    (*H).push_back((*pts)[i]);
  }
  int L = (*H).size();
  (*H).push_back((*pts)[(*pts).size()-2]);
  for(int i=(*pts).size()-3; i>=0; i=i-1) {
    while ((*H).size()>=L+1 && cross(fromto((*H)[(*H).size()-2],(*H)[(*H).size()-1]),
				  fromto((*H)[(*H).size()-1],(*pts)[i]))>=0)
      (*H).pop_back();
    (*H).push_back((*pts)[i]);
  }
  (*H).pop_back();
  return H;
}


int main(int argc, char **argv) {
  
  if (argc==1) {
    printf("This program should be invoked with number of points given on startup!\n");
    printf("For example, if the executable name is 'elliptiser' and you want to elliptise 14 points, enter:\n");
    printf("$ ./elliptiser 14\n");
    return 0;
  }

  if (strcmp(argv[1],"--help")==0) {
    printf("Finds smallest-area ellipse covering given points\n");
    printf("Usage:\n");
    printf("$ elliptiser N\n");
    printf("where N is the number of points you're planning to click\n");
    return 0;
  }
  
  vector<point> *given = new vector<point>;
  bool collected = false;
  int N; sscanf(argv[1],"%d",&N);
  vector<point> *filtrd;


  sf::ContextSettings settings;
  settings.antialiasingLevel = 8;
  sf::RenderWindow window(sf::VideoMode(1200, 800), "My Little Elliptiser",sf::Style::Default, settings);

  
  sf::RectangleShape Ox(sf::Vector2f(1150.f,2.f));
  Ox.setFillColor(sf::Color::Black);
  Ox.setPosition(20.0,399.0);
  sf::ConvexShape tipX;
  tipX.setPointCount(3);
  tipX.setPoint(0,sf::Vector2f(1162.f,403.f));
  tipX.setPoint(1,sf::Vector2f(1177.f,400.f));
  tipX.setPoint(2,sf::Vector2f(1162.f,397.f));
  tipX.setFillColor(sf::Color::Black);
  
  sf::RectangleShape Oy(sf::Vector2f(2.f,750.f));
  Oy.setFillColor(sf::Color::Black);
  Oy.setPosition(599.0,20.0);
  sf::ConvexShape tipY;
  tipY.setPointCount(3);
  tipY.setPoint(0,sf::Vector2f(603.f,27.f));
  tipY.setPoint(1,sf::Vector2f(600.f,12.f));
  tipY.setPoint(2,sf::Vector2f(597.f,27.f));
  tipY.setFillColor(sf::Color::Black);


  sf::VertexArray gridx(sf::Lines,22);
  for (int i=-5; i<=5; i=i+1) {
    if (i!=0) {
      int j=i+5;
      gridx[2*j].position = {600.f+100.f*i,20.f};
      gridx[2*j].color = sf::Color{127,127,127,255};
      gridx[2*j+1].position = {600.f+100.f*i,770.f};
      gridx[2*j+1].color = sf::Color{127,127,127,255};
      }
  }
  sf::VertexArray gridy(sf::Lines,18);
  for (int i=-3; i<=3; i=i+1) {
    if (i!=0) {
      int j=i+3;
      gridy[2*j].position = {20.f,400.f+100.f*i};
      gridy[2*j].color = sf::Color{127,127,127,255};
      gridy[2*j+1].position = {1170.f,400.f+100.f*i};
      gridy[2*j+1].color = sf::Color{127,127,127,255};
      }
  }


  vector<sf::CircleShape> ptz;
  vector<sf::CircleShape> hptz;
  sf::VertexArray *llps;
  
  while (window.isOpen())
  {
      sf::Event event;
      while (window.pollEvent(event))
      {
	if (event.type == sf::Event::Closed) {
	  if (collected)
	    delete filtrd;
	   window.close();
	}
	 if (!collected)
	   if (event.type == sf::Event::MouseButtonPressed)
	     if (event.mouseButton.button == sf::Mouse::Left) {
	       float x,y;
	       x=event.mouseButton.x; y=event.mouseButton.y;
	       sf::CircleShape pd(3.f);
	       pd.setFillColor(sf::Color::Red);
	       pd.setOrigin(3.f,3.f);
	       pd.setPosition(x,y);
	       ptz.push_back(pd);
	       point pm;
	       pm.x=x; pm.y=-y;
	       given->push_back(pm);
	       if (given->size()==N) {
		 collected=true;
		 filtrd = hull(given);
		 delete given;
		 for (int i=0; i<filtrd->size(); i=i+1) {
		   sf::CircleShape hpd(4.f);
		   hpd.setFillColor(sf::Color::Green);
		   hpd.setOrigin(4.f,4.f);
		   hpd.setPosition((*filtrd)[i].x,-(*filtrd)[i].y);
		   hptz.push_back(hpd);
		 }
		 llps = drellipse(encloser(filtrd),150);
	       }
	       
	   }
      }

      window.clear(sf::Color::White);
      window.draw(Ox);
      window.draw(tipX);
      window.draw(Oy);
      window.draw(tipY);
      window.draw(gridx);
      window.draw(gridy);
      for (int i=0; i<ptz.size(); i=i+1)
	window.draw(ptz[i]);
      if (collected)
	for (int i=0; i<hptz.size(); i=i+1)
	  window.draw(hptz[i]);
      if (collected)
	window.draw(*llps);
      window.display();
    }

  return 0;
}
