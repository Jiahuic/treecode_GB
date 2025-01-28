
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "global.h"
#include "tree.h"

/*
 * global variables
 */
double Ra=1.0, l2Inv[3]; /* axes of sphere */
panel *mkIco(int lev, int *nPanels, ssystem *sys);
void liftOntoEllipsoid(panel *snglist, double a, double b, double c, int nSurf);

int main(int nargs, char *argv[]) {
  int nPnls, refineLev = 3, order=-1, orderMom=0;
  ssystem *sys;
  panel *inputLst;

  CALLOC(sys, 1, ssystem, ON, AMISC);
  sys->height = 2;
  sys->depth = -1;
  sys->maxSepRatio = 0.8;
  sys->maxQuadOrder = 1;
  sys->depth = 4;

  for ( int i = 1; i < nargs; i++ ) {
    if ( argv[i][0] == '-' ) {
      switch ( argv[i][1] ) {
        case 'R':
          sscanf(argv[i], "-R%lf", &Ra);
          break;
        case 'L':
          sscanf(argv[i], "-L%d", &refineLev);
          break;
      }
    }
  }
  printf("Ra = %f\n", Ra);
  printf("refineLev = %d\n", refineLev);

  if ( sys->depth < 0 ) {
    printf("Select tree depth > ");
    scanf("%d", &sys->depth);
    if( sys->depth < 1 ) {
      printf("bad tree depth: %d\n", sys->depth );
      exit(0);
    }
  }

  print_ssystem(sys);

  l2Inv[0] = 1.0/SQR(Ra);
  l2Inv[1] = 1.0/SQR(Ra);
  l2Inv[2] = 1.0/SQR(Ra);

  inputLst = mkIco(refineLev, &nPnls, sys);
  liftOntoEllipsoid(inputLst, Ra, Ra, Ra, 1);

  printf("--- nPnls=%d nLev=%d ord=%d SepRat=%lg qOrd=%d\n",
         nPnls, sys->depth, order, sys->maxSepRatio, sys->maxQuadOrder );

  sys->pnlOLst = inputLst;
  tree_init(sys, inputLst, order, orderMom);

  cube *cb;
  dispCube(sys->cubeList[0]);
  for (cb=sys->cubeList[2]; cb!=NULL; cb=cb->next){
    dispCube(cb);
  }
}
