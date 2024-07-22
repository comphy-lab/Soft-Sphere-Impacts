

#include "axi.h"
#include "navier-stokes/centered.h"
#include "fractions.h"
#include "tag.h"
#include "heights.h"

scalar f[];

char filename[80], name[80];

int main(int a, char const *arguments[]) {
  sprintf (filename, "%s", arguments[1]);
  sprintf(name, "%s", arguments[2]);
  restore (file = filename);
  // fprintf(ferr, "Ohd %3.2e, We %g\n", Ohd, We);
  // return 1;

  // boundary conditions
  u.t[left] = dirichlet(0.0);
  f[left] = dirichlet(0.0);

  f.prolongation = refine_bilinear;
  boundary((scalar *){f, u.x, u.y, p});

  // tag all liquid parts starts
  scalar d[];
  double threshold = 1e-4;
  foreach(){
    d[] = (f[] > threshold);
  }
  int n = tag (d), size[n];
  for (int i = 0; i < n; i++){
    size[i] = 0;
  }
  foreach_leaf(){
    if (d[] > 0){
      size[((int) d[]) - 1]++;
    }
  }
  int MaxSize = 0;
  int MainPhase = 0;
  for (int i = 0; i < n; i++){
    // fprintf(ferr, "%d %d\n",i, size[i]);
    if (size[i] > MaxSize){
      MaxSize = size[i];
      MainPhase = i+1;
    }
  }
  // tag all liquid parts ends

  scalar sf[];
  foreach()
    sf[] = (4.*f[] + 
	    2.*(f[0,1] + f[0,-1] + f[1,0] + f[-1,0]) +
	    f[-1,-1] + f[1,-1] + f[1,1] + f[-1,1])/16.;
  sf.prolongation = refine_bilinear;
  boundary ({sf});
  /*
  Do calculations start
  */
  face vector s[];
  s.x.i = -1;
  double yMax = 0.0;
  
  foreach(){
    if (f[] > 1e-6 && f[] < 1. - 1e-6 && d[] == MainPhase) {
      coord n1 = facet_normal (point, f, s);
      double alpha1 = plane_alpha (f[], n1);
      coord segment1[2];
      if (facets (n1, alpha1, segment1) == 2){
        double x1 = x + (segment1[0].x+segment1[1].x)*Delta/2.;
        double y1 = y + (segment1[0].y+segment1[1].y)*Delta/2.;
        if (x1>0 && x1<0.01){
            if (y1 > yMax){
                 yMax = y1;
            }
        }
      }
    }
  }

  /*
  Do calculations end
  */
//Write file
  FILE *fp;
  fp = fopen (name, "a");
  restore (file = filename);

  if (t == 0){
    fprintf(ferr, "t Rfoot\n");
    fprintf(fp, "t Rfoot\n");    
  }

  fprintf(ferr, "%6.5e %6.5e\n", t, yMax);
  fprintf(fp, "%6.5e %6.5e\n", t, yMax);
  fclose(fp);
}
