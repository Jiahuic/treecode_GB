#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "global.h"
#include "tree.h"

int nSurf = 0;

/*
 * Uniform refinement: Each panel is refined into four panels of equal size.
 * The refined panels are inserted into the panel list.
 * Only triangles are supported
 */
void unifRefine(panel *pnlList, int *nPanels, double raduis) {
   panel *pnl, *pnlNext;
   int j, k, nSurf;
   double v0[3], v1[3], v2[3], w0[3], w1[3], w2[3];
   double leng1, leng2, leng3;

   pnl=pnlList;

   while ( pnl != NULL ) {

     pnlNext = pnl->next;
     if ( pnl->shape != 3 ) {
       printf("\n**error** unifRefine(): only triangles\n");
       exit(1);
     }

     leng1 = leng2 = leng3 = 0.0;
     for ( k=0; k<3; k++ ) {
       v0[k] = pnl->vtx[0][k];
       v1[k] = pnl->vtx[1][k];
       v2[k] = pnl->vtx[2][k];
       w0[k] = 0.5*(v1[k] + v2[k]);
       w1[k] = 0.5*(v2[k] + v0[k]);
       w2[k] = 0.5*(v0[k] + v1[k]);
       leng1 += SQR(w0[k]);
       leng2 += SQR(w0[k]);
       leng3 += SQR(w0[k]);
     }
     leng1 = sqrt(leng1);
     leng2 = sqrt(leng2);
     leng3 = sqrt(leng3);

     for ( k=0; k<3; k++ ) {
       w0[k] *= raduis/leng1;
       w1[k] *= raduis/leng2;
       w2[k] *= raduis/leng3;
       pnl->vtx[0][k] = v0[k];
       pnl->vtx[1][k] = w2[k];
       pnl->vtx[2][k] = w1[k];
     }

     CALLOC(pnl->next, 1, panel, ON, APVE);
     pnl = pnl->next;
     for ( k=0; k<3; k++ ) {
       pnl->vtx[0][k] = w2[k];
       pnl->vtx[1][k] = v1[k];
       pnl->vtx[2][k] = w0[k];
     }
     pnl->nSurf = nSurf++;
     //printf("%d\n",pnl->nSurf);
     pnl->shape = 3;

     CALLOC(pnl->next, 1, panel, ON, APVE);
     pnl = pnl->next;
     for ( k=0; k<3; k++ ) {
       pnl->vtx[0][k] = w0[k];
       pnl->vtx[1][k] = v2[k];
       pnl->vtx[2][k] = w1[k];
     }
     pnl->nSurf = nSurf++;
     //printf("%d\n",pnl->nSurf);
     pnl->shape = 3;

     CALLOC(pnl->next, 1, panel, ON, APVE);
     pnl = pnl->next;
     for ( k=0; k<3; k++ ) {
       pnl->vtx[0][k] = w0[k];
       pnl->vtx[1][k] = w1[k];
       pnl->vtx[2][k] = w2[k];
     }
     pnl->nSurf = nSurf++;
     //printf("%d\n",pnl->nSurf);
     pnl->shape = 3;

     pnl->next = pnlNext;
     pnl = pnlNext;
   }

   *nPanels *= 4;

}/* unifRefine */

/*
 * Make the icosahedron.
 * The panels are oriented such that the normal determined by the right hand
 * rule points into the exterior of the cube.
 *
 * Parameters:
 *  lev       refinement level
 *  nPanels   (output) number of panels:
 *                     lev=0=>nPanels=48, lev=1=>nPanels=192,...
 */
panel *mkIco(int lev, int *nPanels, ssystem *sys) {
  panel *pnlList=NULL, *pnl;

  int cnt, j, k, level;
  /* 20 panels of initial triangulation */
  double vtx[20][9] = {
    {0.85065080835203993218, 0, 0.525731112119133606025,
      0.85065080835203993218, 0, -0.525731112119133606025,
      0.525731112119133606025, 0.85065080835203993218},
    {0.85065080835203993218, 0, 0.525731112119133606025,
      0.525731112119133606025, -0.85065080835203993218, 0,
      0.85065080835203993218, 0, -0.525731112119133606025},
    {-0.85065080835203993218, 0, -0.525731112119133606025,
      -0.85065080835203993218, 0, 0.525731112119133606025,
      -0.525731112119133606025, 0.85065080835203993218},
    {-0.85065080835203993218, 0, 0.525731112119133606025,
      -0.85065080835203993218, 0, -0.525731112119133606025,
      -0.525731112119133606025, -0.85065080835203993218},
    {0.525731112119133606025, 0.85065080835203993218, 0,
      -0.525731112119133606025, 0.85065080835203993218, 0,
      0, 0.525731112119133606025, 0.85065080835203993218},
    {-0.525731112119133606025, 0.85065080835203993218, 0,
      0.525731112119133606025, 0.85065080835203993218, 0,
      0, 0.525731112119133606025, -0.85065080835203993218 },
    {0, -0.525731112119133606025, -0.85065080835203993218,
      0, 0.525731112119133606025, -0.85065080835203993218,
      0.85065080835203993218, 0, -0.525731112119133606025 },
    {0, 0.525731112119133606025, -0.85065080835203993218,
      0, -0.525731112119133606025, -0.85065080835203993218,
      -0.85065080835203993218, 0, -0.525731112119133606025},
    {0.525731112119133606025, -0.85065080835203993218, 0,
      -0.525731112119133606025, -0.85065080835203993218, 0,
      0, -0.525731112119133606025, -0.85065080835203993218},
    {-0.525731112119133606025, -0.85065080835203993218, 0,
      0.525731112119133606025, -0.85065080835203993218, 0,
      0, -0.525731112119133606025, 0.85065080835203993218},
    {0, 0.525731112119133606025, 0.85065080835203993218, 0,
      -0.525731112119133606025, 0.85065080835203993218,
      0.85065080835203993218, 0, 0.525731112119133606025 },
    {0, -0.525731112119133606025, 0.85065080835203993218,
      0, 0.525731112119133606025, 0.85065080835203993218,
      -0.85065080835203993218, 0, 0.525731112119133606025},
    {0.525731112119133606025, 0.85065080835203993218, 0,
      0.85065080835203993218, 0, -0.525731112119133606025,
      0, 0.525731112119133606025, -0.85065080835203993218},
    {0.85065080835203993218, 0, 0.525731112119133606025,
      0.525731112119133606025, 0.85065080835203993218, 0,
      0, 0.525731112119133606025, 0.85065080835203993218 },
    {-0.85065080835203993218, 0, -0.525731112119133606025,
      -0.525731112119133606025, 0.85065080835203993218, 0,
      0, 0.525731112119133606025, -0.85065080835203993218 },
    {-0.525731112119133606025, 0.85065080835203993218, 0,
      -0.85065080835203993218, 0, 0.525731112119133606025,
      0, 0.525731112119133606025, 0.85065080835203993218 },
    {0.85065080835203993218, 0, -0.525731112119133606025,
      0.525731112119133606025, -0.85065080835203993218, 0,
      0, -0.525731112119133606025, -0.85065080835203993218},
    {0.525731112119133606025, -0.85065080835203993218, 0,
      0.85065080835203993218, 0, 0.525731112119133606025,
      0, -0.525731112119133606025, 0.85065080835203993218 },
    {-0.85065080835203993218, 0, -0.525731112119133606025,
      0, -0.525731112119133606025, -0.85065080835203993218,
      -0.525731112119133606025, -0.85065080835203993218, 0},
    {-0.85065080835203993218, 0, 0.525731112119133606025,
      -0.525731112119133606025, -0.85065080835203993218, 0,
      0, -0.525731112119133606025, 0.85065080835203993218}
  };

  /* The unit charge is located at the center of the sphere with some raduis */
  sys->nChar = 1;
  CALLOC(sys->pos, 3*sys->nChar, double, ON, APVE);
  CALLOC(sys->chr, sys->nChar, double, ON, APVE);
  sys->pos[0] = 0.0;
  sys->pos[1] = 0.0;
  sys->pos[2] = 0.0;
  sys->chr[0] = 1.0;
  *nPanels = 20;
  nSurf = 0;

  for ( cnt=0; cnt<20; cnt++ ) {
    /* Allocate panel struct to fill in. */
    if(pnlList == NULL) {
      CALLOC(pnlList, 1, panel, ON, APVE);
      pnl = pnlList;
    }
    else {
      CALLOC(pnl->next, 1, panel, ON, APVE);
      pnl = pnl->next;
    }

    /* Fill in corners. */
    for(j=0;j<3;j++) {
      for(k=0;k<3;k++) {
        pnl->vtx[j][k] = vtx[cnt][j*3+k];
      }
    }
    pnl->shape = 3;
    pnl->nSurf = nSurf++;
    //printf("mkIco%d\n",pnl->nSurf);
  }

  for ( level=1; level<=lev; level++ ){
    unifRefine(pnlList, nPanels, 1.0);
  }

  return pnlList;
} /* mkIco */

/*
 * maps input panles on the ellipsoid
 *      (x/a)^2 + (y/b)^2 + (z/c)^2 = 1
 *
 */
void liftOntoEllipsoid(panel *pnlList, double a, double b, double c, int nSurf) {
  panel *pnl;
  double len;
  int i;

  for ( pnl=pnlList; pnl != NULL; pnl=pnl->next ) {
    for ( i=0; i<pnl->shape; i++ ) {
      /* this gives a better projection */
      len = SQR(pnl->vtx[i][0]) + SQR(pnl->vtx[i][1]) + SQR(pnl->vtx[i][2]);
      len = sqrt(len);
      pnl->vtx[i][0] *= a/len;
      pnl->vtx[i][1] *= b/len;
      pnl->vtx[i][2] *= c/len;
    }
  }
} /* liftOntoEllipsoid */
