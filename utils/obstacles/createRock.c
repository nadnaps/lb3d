/* 
=========================================================================

 Copyright 1999-2012, Owners retain copyrights to their respective works.

 This file is part of lb3d.

 lb3d is free software: you can redistribute it and/or modify it under
 the terms of the GNU Lesser General Public License as published by
 the Free Software Foundation, either version 3 of the License, or (at
 your option) any later version.

 lb3d is distributed in the hope that it will be useful, but WITHOUT
 ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
 FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public
 License for more details.

 You should have received a copy of the GNU Lesser General Public
 License along with lb3d. If not, see <http://www.gnu.org/licenses/>.

=========================================================================

createRock tool 2 v.0.0 
date: 06.07.10
Christian Kunert
Ariel Narváez

The createRock tool generates a xdr domain

To generate empty space use rock color 0.0
To generate solid space use rock color 5.0

the input data has the format:
1) In the beginning a header is needed that contains the dimension in x,y, and z 
     8        # x-dimension
     128      # y-dimension
     512      # z-dimension


2) The color of the domain just to begin with the construction of the structure
     0.0  # domain color at the beginning
# sometimes is usefull to start being the whole domain a rock (5.00)
# and sometimes being the whole domain fluid space (0.00) or even other color. 


3) Domain construction, setting the geometries (wall, cube_trj, rctgl_trj, sphr_trj, 
cube, sphere, ellipse, sqr_column, or circ_column)
composing the domain. They can be set at any order and as many times is desired.


If a wall is being set, this crosses the whole domain and its parameters are: 
     1 3 1    # starting point (int)
     0 4 0    # wall growing direction (int)
#  in the 'wall growing direction'
#  can only be a vector with one component different from zero, 
#  this no-zero component is the direction of the wall's width.
#  a positive value means that the wall width grows in the direction of the depth,
#  negative the wall width grows in the other direction
     0.0      # colour (float)


If a cube_trj is being set, this is contructucted by a cube oriented in the x, y, 
and z axes, this cube travels from the starting to the ending point. 
The center of the cube moves along the trajectory.
The parameters are:
     7 67 0   # starting point (int) 
     8 61 512 # ending point (int)
     66       # width, centered in the line that connect 
              # the starting and ending points 
     0.0      # colour (float)


If a rctgl_trj is being set, this is contructucted by a rectangle oriented in the x, y, 
and z axes, this rectangle travels from the starting to the ending point. 
The center of the rectangle moves along the trajectory.
The parameters are:
     7 67 0   # starting point (int) 
     8 61 512 # ending point (int)
     66 10 20 # dimensions of the rectangle, centered in the line that connect 
              # the starting and ending points 
     0.0      # colour (float)


If a sphr_trj is being set, this is contructucted by a sphere of radious, 
this sphere travels from the starting to the ending point. 
The center of the sphere moves along the trajectory.
The parameters are:
     7 67 0   # starting point (int) 
     8 61 512 # ending point (int)
     6.6      # radious, centered in the line that connect 
              # the starting and ending points (float)
     0.0      # colour (float)


If a cube is being set, its parameters are:
     5 6 8   # center (int)
     15      # width (int)
     0.0     # colour (float)


If a sphere is being set, its parameters are:
     4.5 5.1 10.10   # center point (float)
     3.2             # sphere radious (float)
     0.0             # colour (float)


If ellipse, its parameters are:
     5.1 7.1 10.9   # focus point (float)
     8.5 5.3 10.4   # focus point (float)
     4.2            # sum of the distances to the two focus points (float)
     0.0            # colour (float)


If a sqr_column is being set, a column with square cross_section is being build,
this can only be oriented in the direction of the axis x, y, and z.
Its parameters are:
     5 7 10   # starting center point (int) 
     0 0 10   # direction an length (int)
#  this last parameter can only be a vector with one component different from zero.
     4.2      # sqr_column width (int)
     0.0      # colour (float)


If a circ_column is being set, a column with circular cross_section is being build, 
this column can only be oriented in the direction of the axis x, y, and z.
Its parameters are:
     5.2 7.2 10.2   # starting center point (float)
     0 0 10         # direction an length (int)
#  this last parameter can only be a vector with one component different from zero.
     4.2            # circ_column radious (float)
     0.0            # colour (float)


4) Write 'end' to stop reading the geometries, from this time one the output files will be written.



compile the program:

gcc -lm createRock.c -o createRock

running:
./createRock input-data

output: input-data.txt input-data.xtr

*/

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <rpc/rpc.h>
#include <rpc/xdr.h>
#include <math.h>
int f_domain(int c1, int cn);
int f_trj(int c1, int c2, int pmax, int p);
float f_sphere(float cx, float cy, float cz, int i, int j, int k);
float f_circ(float cx, float cy, int i, int j);
float f_ellipse(float cx1, float cy1, float cz1, float cx2, float cy2, float cz2, int i, int j, int k, float radious);

#define MAXSTRING 160

int main (int argc, char *argv[]){
  FILE	    *f0, *f1, *f2;
  int       i, j, k;
  XDR	    xdrs;
  char      *fname;
  char      extxdr[] = ".xdr", ext[] = ".txt";
  char      xdrfile[MAXSTRING], txtfile[MAXSTRING];
  char      command[MAXSTRING];
  int       nx, ny, nz;
  float     ***pos;

  char      shape[] = "begin";

  float	    colour;

  int       wll; /* wall number */
  int       wi, wj, wk; /* wall position */
  int       vg_i, vg_j, vg_k;  /* wall growing direction */

  int       cb_trj; /* cube_trj number */
  int       ti1, ti2, tj1, tj2, tk1, tk2 ; /* cube_trj trajectory */
  int       width;  /* width */
  int       widthr, widthl, ii, jj, kk, t, tmax;

  int       rctgl_trj; /* rctgl_trj number */
  int       width_x, width_y, width_z;  /* width */
  int       width_xr, width_xl, width_yr, width_yl, width_zr, width_zl;

  int       sp_trj; /* sphr_trj number */
  float     radious;  /* radious */

  int       cb; /* cube number */

  int       sphr; /* sphere number */
  float     sx, sy, sz; /* center sphere */

  int       elps; /* ellipse number */
  float     ellx1, ellx2, elly1, elly2, ellz1, ellz2 ; /* ellipse */
  int       ti_max, ti_min, tj_max, tj_min, tk_max, tk_min;

  int       sqr_clm; /* sqr_column number */

  int       circ_clm; /* circ_column number */

  if (argc!=2){
    fprintf(stderr,"Syntax: %s <filename> \n",argv[0]);
    return -1;
  }
  
  fname=argv[1];
  
  /* Input file */
  f0 = fopen(fname,"r");
  /* 3D grid */
  fscanf(f0,"%d",&nx);
  fscanf(f0,"%d",&ny);
  fscanf(f0,"%d",&nz);  
  printf("domain dimension %d %d %d \n", nx, ny, nz);


  /* Reading color domain */
  fscanf(f0,"%f",&colour); 
  printf("starting color domain %f \n", colour);


  sprintf(txtfile,"%s%s",fname,ext);

  /* Allocate tensor (i.e. 3d array) and set pointers */
  pos=(float ***)malloc(nx*sizeof(float **));
  if(!pos){
    fprintf(stderr,"error: (1) can't allocate memory \n");
    exit(1);
  }
  pos[0]=(float **)malloc(nx*ny*sizeof(float *));
  if(!pos[0]){
    fprintf(stderr,"error: (2) can't allocate memory \n"); 
    exit(1);
  }
  pos[0][0]=(float *)calloc(nx*ny*nz,sizeof(float));
  if(!pos[0][0]){
    fprintf(stderr,"error: (3) can't allocate memory \n");
    exit(1);
  }
  
  for( j=1; j<ny; j++ ) {
    pos[0][j]=pos[0][j-1]+nz;
  }

  for( i=1; i<nx; i++ ) {
    pos[i]=pos[i-1]+ny;
    pos[i][0]=pos[i-1][0]+ny*nz;
    for( j=1; j<ny; j++ ) {
      pos[i][j]=pos[i][j-1]+nz;
    }
  } 
   
  sprintf(txtfile,"%s%s",fname,ext);
  sprintf(xdrfile,"%s%s",fname,extxdr);
  /* Open ASCII file */
  f1 = fopen(txtfile,"w");
  /* Create xdrfile, stop if it already exists/overwrite? */
  sprintf(command,"touch %s ",xdrfile);
  system(command);
  if(NULL==(f2=fopen(xdrfile,"w"))){
    perror("unable to open output file");
    return -1;
  }        
  printf("%s \n",xdrfile);

  /* In the beginning set all points to color_dmn */ 
  for ( k=0; k<nz; ++k ) {
    for ( j=0; j<ny; ++j ) {
      for ( i=0; i<nx; ++i ) {
	pos[i][j][k] = colour;
      }
    }
  }

  /* Reading geometry parameters */
  
  printf("start reading and creating the geometries\n");

  wll = 0;
  cb_trj = 0;
  rctgl_trj = 0;
  sp_trj = 0;
  cb = 0;
  sphr = 0;
  elps = 0;
  sqr_clm = 0;
  circ_clm = 0;

  while ( strcmp(shape,"end") != 0 ) {
    fscanf(f0,"%s",shape);

    if ( strcmp(shape,"wall") == 0 ) {
      fscanf(f0, "%d %d %d", &wi, &wj, &wk);
      fscanf(f0, "%d %d %d", &vg_i, &vg_j, &vg_k);
      fscanf(f0, "%f", &colour);
      wll += 1;   
      printf("wall %d: %d %d %d %d %d %d %f\n", wll, wi, wj, wk, vg_i, vg_j, vg_k, colour);
     
      /* set the wall */     
      if ( vg_i!=0 && vg_j==0 && vg_k==0 ) {
	if ( vg_i > 0 ) {
	  for ( k=0; k<nz; ++k ) {
	    for ( j=0; j<ny; ++j ) {
	      for ( i=wi; i<=wi+vg_i-1; ++i ) {
		ii = f_domain(i,nx);
		pos[ii][j][k] = colour; 
	      }
	    }
	  }
	}
	else {
	  for ( k=0; k<nz; ++k ) {
	    for (j=0; j<ny; ++j ) {
	      for (i=wi+vg_i+1; i<=wi; ++i ) {
		ii = f_domain(i,nx);
		pos[ii][j][k] = colour; 
	      }
	    }
	  }
	}
      }
      else if ( vg_i==0 && vg_j!=0 && vg_k==0 ) {
	if ( vg_j>0 ) {
	  for ( k=0; k<nz; ++k ) {
	    for (j=wj; j<=wj+vg_j-1 ;++j ) {
	      for (i=0 ; i<nx; ++i) {
		jj = f_domain(j,ny);
		pos[i][jj][k] = colour; 
	      }
	    }
	  }
	}
	else {
	  for ( k=0; k<nz; ++k ) {
	    for ( j=wj+vg_j+1; j<=wj; ++j ) {
	      for ( i=0; i<nx; ++i ) {
		jj = f_domain(j,ny);
		pos[i][jj][k] = colour; 
	      }
	    }
	  }
	}
      }  
      else if ( vg_i==0 && vg_j==0 && vg_k!=0 ) {
	if ( vg_k>0 ) {
	  for ( k=wk; k<=wk+vg_k-1; ++k ) {
	    for ( j=0; j<ny; ++j) {
	      for ( i=0; i<nx; ++i) {
		kk = f_domain(k,nz);
		pos[i][j][kk] = colour; 
	      }
	    }
	  }
	}
	else {
	  for ( k=wk+vg_k+1; k<=wk; ++k ) {
	    for ( j=0; j<ny; ++j ) {
	      for ( i=0; i<nx; ++i ) {
		kk = f_domain(k,nz);
		pos[i][j][kk] = colour; 
	      }
	    }
	  }
	}
      }
      else {
	printf("Attention, wrong definition of the wall %d\n",wll);
	exit(1);
      }   
      /* end of wall */

    } else if ( strcmp(shape,"cube_trj") == 0 ) { 
      fscanf(f0, "%d %d %d", &ti1, &tj1, &tk1);
      fscanf(f0, "%d %d %d", &ti2, &tj2, &tk2);
      fscanf(f0, "%d", &width);
      fscanf(f0, "%f", &colour);
      cb_trj += 1; 
      printf("cube_trj %d: %d %d %d %d %d %d %d %f\n", cb_trj, ti1, tj1, tk1, ti2, tj2, tk2, width, colour);

      /* cube_trj with 1 node width */ 
      if ( width >= 1 ) {

	/* t is the parameter to contruc the line between (x1,y1,z1) and (x2,y2,z2) */
	tmax = abs(ti1-ti2); 
	if ( abs(tj1-tj2)>=tmax ) {
	  tmax = abs(tj1-tj2);
	}
	if ( abs(tk1-tk2)>=tmax ) {
	  tmax = abs(tk1-tk2);
	}
	
	if ( width == 1 ) {
	  for ( t=0; t<=tmax; ++t) {
	    /* tmax = 0 */
	    if ( tmax == 0 ) {
	      i = ti1;
	      j = tj1;	
	      k = tk1;
	    }
	    else {
	      i = f_trj(ti1, ti2, tmax, t);
	      j = f_trj(tj1, tj2, tmax, t);	
	      k = f_trj(tk1, tk2, tmax, t);
	    }
	    pos[i][j][k] = colour; 
	  }
	}
	else {
	  /* width is odd */
	  if ( width/2 != (width+1)/2 ) {
	    widthr = (width-1)/2;
	    widthl = (width-1)/2; 
	  }
	  /* width is even */
	  else {
	    widthr = width/2;
	    widthl = width/2 - 1;
	  }

	  printf("%d %d %d\n", width, widthr, widthl);
	  for ( t=0; t<=tmax; ++t ) {
	    /* tmax = 0 */
	    if (tmax == 0){
	      i = ti1;
	      j = tj1;	
	      k = tk1;
	    }
	    else {
	      i = f_trj(ti1,ti2,tmax,t);
	      j = f_trj(tj1,tj2,tmax,t);	
	      k = f_trj(tk1,tk2,tmax,t);
	    }
	    for ( ii=-widthl; ii<=widthr; ++ii ) {
	      for ( jj=-widthl; jj<=widthr; ++jj ) {
		for ( kk=-widthl; kk<=widthr; ++kk ) {
		  if (i+ii>=0 && j+jj>=0 && k+kk>= 0 && i+ii<=nx-1 && j+jj<=ny-1 && k+kk<=nz-1) {
		  pos[i+ii][j+jj][k+kk] = colour; 
		  }
		}
	      }
	    }
	  }
	}
      }
      else {
	printf("Attention, wrong definition of the cube_trj %d\n",cb_trj);
	exit(1);
      } 
      /* end of cube_trj */


    } else if ( strcmp(shape,"rctgl_trj") == 0 ) { 
      fscanf(f0, "%d %d %d", &ti1, &tj1, &tk1);
      fscanf(f0, "%d %d %d", &ti2, &tj2, &tk2);
      fscanf(f0, "%d %d %d", &width_x, &width_y, &width_z);
      fscanf(f0, "%f", &colour);
      rctgl_trj += 1; 
      printf("cube_trj %d: %d %d %d %d %d %d %d %d %d %f\n", cb_trj, ti1, tj1, tk1, ti2, tj2, tk2, width_x, width_y, width_z, colour);

      /* rctgl_trj with 1 node width */ 
      if ( width_x >= 1 && width_y >= 1 && width_z >= 1) {

	/* t is the parameter to contruc the line between (x1,y1,z1) and (x2,y2,z2) */
	tmax = abs(ti1-ti2); 
	if ( abs(tj1-tj2)>=tmax ) {
	  tmax = abs(tj1-tj2);
	}
	if ( abs(tk1-tk2)>=tmax ) {
	  tmax = abs(tk1-tk2);
	}
	
	/* width_x is odd */
	if ( width_x/2 != (width_x+1)/2 ) {
	  width_xr = (width_x-1)/2;
	  width_xl = (width_x-1)/2; 
	}
	/* width_x is even */
	else {
	  width_xr = width_x/2;
	  width_xl = width_x/2 - 1;
	}
	
	/* width_y is odd */
	if ( width_y/2 != (width_y+1)/2 ) {
	  width_yr = (width_y-1)/2;
	  width_yl = (width_y-1)/2; 
	  }
	/* width_y is even */
	else {
	  width_yr = width_y/2;
	  width_yl = width_y/2 - 1;
	}
	
	/* width_z is odd */
	if ( width_z/2 != (width_z+1)/2 ) {
	  width_zr = (width_z-1)/2;
	  width_zl = (width_z-1)/2; 
	}
	/* width_z is even */
	else {
	  width_zr = width_z/2;
	  width_zl = width_z/2 - 1;
	}

	for ( t=0; t<=tmax; ++t ) {
	  /* tmax = 0 */
	  if (tmax == 0){
	    i = ti1;
	    j = tj1;	
	    k = tk1;
	    }
	  else {
	    i = f_trj(ti1,ti2,tmax,t);
	    j = f_trj(tj1,tj2,tmax,t);	
	    k = f_trj(tk1,tk2,tmax,t);
	  }
	  for ( ii=-width_xl; ii<=width_xr; ++ii ) {
	    for ( jj=-width_yl; jj<=width_yr; ++jj ) {
	      for ( kk=-width_zl; kk<=width_zr; ++kk ) {
		if (i+ii>=0 && j+jj>=0 && k+kk>= 0 && i+ii<=nx-1 && j+jj<=ny-1 && k+kk<=nz-1) {
		  pos[i+ii][j+jj][k+kk] = colour; 
		}
	      }
	    }
	  }
	}
      }
      else {
	printf("Attention, wrong definition of the rctgl_trj %d\n",rctgl_trj);
	exit(1);
      } 
      /* end of rctgl_trj */
	
    } else if ( strcmp(shape,"sphr_trj") == 0 ) { 
      fscanf(f0, "%d %d %d", &ti1, &tj1, &tk1);
      fscanf(f0, "%d %d %d", &ti2, &tj2, &tk2);
      fscanf(f0, "%f", &radious);
      fscanf(f0, "%f", &colour);
      sp_trj += 1; 
      printf("sphr_tr %d: %d %d %d %d %d %d %f %f\n", sp_trj, ti1, tj1, tk1, ti2, tj2, tk2, radious, colour);

     /* sphr_trj with 1 node width */ 
      if ( radious >= 1.0 ) {

	/* t is the parameter to contruc the line between (x1,y1,z1) and (x2,y2,z2) */
	tmax = abs(ti1-ti2); 
	if (abs(tj1-tj2)>=tmax){
	  tmax = abs(tj1-tj2);
	}
	if (abs(tk1-tk2)>=tmax){
	  tmax = abs(tk1-tk2);
	}
	
	for ( t=0; t<=tmax; ++t ) {
	  /* tmax = 0 */
	  if (tmax == 0){
	    i = ti1;
	    j = tj1;	
	    k = tk1;
	  } else {
	    i = f_trj(ti1, ti2, tmax, t);
	    j = f_trj(tj1, tj2, tmax, t);	
	    k = f_trj(tk1, tk2, tmax, t);
	  }
	  width = radious + 1;
	  for ( ii=-width; ii<=width; ++ii ) {
	    for ( jj=-width; jj<=width; ++jj ) {
	      for ( kk=-width; kk<=width; ++kk ) {
		if (i+ii>=0 && j+jj>=0 && k+kk>= 0 && i+ii<=nx-1 && j+jj<=ny-1 && k+kk<=nz-1) {
		  if ( pow(radious,0.2e1) >= f_sphere((float)i+(float)ii, (float)j+(float)jj, (float)k+(float)kk, i, j, k) ) {
		    pos[i+ii][j+jj][k+kk] = colour; 
		  }
		}
	      }
	    }
	  }
	}
      }
      else {
	printf("Attention, wrong definition of the sphr_trj %d\n",sp_trj);
	exit(1);
      } 
      /* end of sphr_trj */
      
    } else if ( strcmp(shape,"cube") == 0 ) {
      fscanf(f0, "%d %d %d", &wi, &wj, &wk);
      fscanf(f0, "%d", &width);
      fscanf(f0, "%f", &colour);
      cb += 1;
      printf("cube %d: %d %d %d %d %f \n", cb, wi, wj, wk, width, colour);

      /* width is odd */
      if ( width/2 != (width+1)/2 ) {
	widthr = (width-1)/2;
	widthl = (width-1)/2; 
      }
      /* width is even */
      else {
	widthr = width/2;
	widthl = width/2 - 1;
      }

      /* set all points to colour_sphe 3d */
      for ( k=wk-widthl; k<=wk+widthr; ++k) {
        for ( j=wj-widthl; j<=wj+widthr; ++j) {
	  for ( i=wi-widthl; i<=wi+widthr; ++i) {
	    ii = f_domain(i, nx);
	    jj = f_domain(j, ny);
	    kk = f_domain(k, nz);
	    pos[ii][jj][kk] = colour;
	  }
        }
      }
      /* end of cube */

    } else if ( strcmp(shape,"sphere") == 0 ) { 
      fscanf(f0, "%f %f %f", &sx, &sy, &sz);
      fscanf(f0, "%f", &radious);
      fscanf(f0, "%f", &colour);
      sphr += 1;
      printf("sphere %d: %f %f %f %f %f \n", sphr, sx, sy, sz, radious, colour);
      
      width = radious +2; 
      wi = sx;
      wj = sy;
      wk = sz;

      /* set all points to colour_sphe 3d */
      for ( k=wk-width; k<=wk+width; ++k) {
	for ( j=wj-width; j<=wj+width; ++j) {
	 for ( i=wi-width; i<=wi+width; ++i) {
	   ii = f_domain(i, nx);
	   jj = f_domain(j, ny);
	   kk = f_domain(k, nz);
	   if( pow(radious,0.2e1) >= f_sphere(sx, sy, sz, ii, jj, kk) ) {
	     pos[ii][jj][kk] = colour;
	   }
	 }
	}
      }
      /* end of spheres */

    } else if ( strcmp(shape,"ellipse") == 0 ) { 
      fscanf(f0, "%f %f %f", &ellx1, &elly1, &ellz1);
      fscanf(f0, "%f %f %f", &ellx2, &elly2, &ellz2);
      fscanf(f0, "%f", &radious);
      fscanf(f0, "%f", &colour);
      elps += 1;
      printf("ellipse %d: %f %f %f %f %f %f %f %f\n", elps, ellx1, elly1, ellz1, ellx2, elly2, ellz2, radious, colour);
 
      if ( ellx1 >= ellx2 ) {
	ti_max = (int)ellx2 + (int)radious;
	ti_max = f_domain(ti_max,nx);
	ti_min = (int)ellx1 - (int)radious;
	ti_min = f_domain(ti_min,nx);
      }
      else {
	ti_max = (int)ellx1 + (int)radious;
	ti_max = f_domain(ti_max,nx);
	ti_min = (int)ellx2 - (int)radious;
	ti_min = f_domain(ti_min,nx);
      }
      
      if ( elly1 >= elly2 ) {
	tj_max = (int)elly2 + (int)radious;
	tj_max = f_domain(tj_max,nx);
	tj_min = (int)elly1 - (int)radious;
	tj_min = f_domain(tj_min,nx);
      }    
      else {
	tj_max = (int)elly1 + (int)radious;
	tj_max = f_domain(tj_max,nx);
	tj_min = (int)elly2 - (int)radious;
	tj_min = f_domain(tj_min,nx);
      }
      
      if ( ellz1 >= ellz2 ) {
	tk_max = (int)ellz2 + (int)radious;
	tk_max = f_domain(tk_max,nx);
	tk_min = (int)ellz1 - (int)radious;
	tk_min = f_domain(tk_min,nx);
      }
      else {
	tk_max = (int)ellz1 + (int)radious;
	tk_max = f_domain(tk_max,nx);
	tk_min = (int)ellz2 - (int)radious;
	tk_min = f_domain(tk_min,nx);
      }
      
      for ( k=tk_min; k<=tk_max; ++k ) {
	for ( j=tj_min; j<=tj_max; ++j ) {
	  for ( i=ti_min; i<=ti_max; ++i ) {
	    if ( radious >= f_ellipse(ellx1, elly1, ellz1, ellx2, elly2, ellz2, i, j, k, radious) ) {
	      pos[i][j][k] = colour; 
	    }
	  }
	}
      }
      /* end of ellipses */

    } else if ( strcmp(shape,"sqr_column") == 0 ) { 
      fscanf(f0, "%d %d %d", &wi, &wj, &wk);
      fscanf(f0, "%d %d %d", &vg_i, &vg_j, &vg_k);
      fscanf(f0, "%d", &width);
      fscanf(f0, "%f", &colour);
      sqr_clm += 1;
      printf("sqr_column %d: %d %d %d %d %d %d %d %f\n", sqr_clm, wi, wj, wk, vg_i, vg_j, vg_k, width, colour); 
      
      if ( width >= 1 ) {
	
	/* width channel is odd */
	if ( width/2 != (width+1)/2 ) {
	  widthr = (width-1)/2;
	  widthl = (width-1)/2; 
	}
	/* width channel is even */
	else{
	  widthr = width/2;
	  widthl = width/2 - 1;
	}
	
	if ( vg_i!=0 && vg_j==0 && vg_k==0 ){
	  if ( vg_i>0 ) {
	    for ( k=wk-widthl; k<=wk+widthr; ++k ) {
	      for ( j=wj-widthl; j<=wj+widthr; ++j ) {
		for ( i=wi; i<=wi+vg_i; ++i ) {
		  ii = f_domain(i, nx);
		  jj = f_domain(j, ny);
		  kk = f_domain(k, nz);
		  pos[ii][jj][kk] = colour; 
		}
	      }
	    }
	  }
	  else {
	    for ( k=wk-widthl; k<=wk+widthr; ++k ) {
	      for ( j=wj-widthl; j<=wj+widthr; ++j ) {
		for ( i=wi+vg_i; i<=wi; ++i ) {
		  ii = f_domain(i, nx);
		  jj = f_domain(j, ny);
		  kk = f_domain(k, nz);
		  pos[ii][jj][kk] = colour; 
		}
	      }
	    }
	  }
	}
	else if ( vg_i==0 && vg_j!=0 && vg_k==0 ){
	  if ( vg_j>0 ) {
	    for ( k=wk-widthl; k<=wk+widthr; ++k ) {
	      for ( j=wj; j<=wj+vg_j; ++j ) {
		for ( i=wi-widthl; i<=wi+widthr; ++i ) {
		  ii = f_domain(i, nx);
		  jj = f_domain(j, ny);
		  kk = f_domain(k, nz);
		  pos[ii][jj][kk] = colour; 
		}
	      }
	    }
	  }
	  else {
	    for ( k=wk-widthl; k<=wk+widthr; ++k ) {
	      for ( j=wj+vg_j; j<=wj; ++j ) {
		for ( i=wi-widthl; i<=wi+widthr; ++i ) {
		  ii = f_domain(i, nx);
		  jj = f_domain(j, ny);
		  kk = f_domain(k, nz);
		  pos[ii][jj][kk] = colour; 
		}
	      }
	    }
	  }
	}
	else if ( vg_i==0 && vg_j==0 && vg_k!=0 ){
	  if ( vg_k>0 ) {
	    for ( k=wk; k<=wk+vg_k; ++k ) {
	      for ( j=wj-widthl; j<=wj+widthr; ++j ) {
		for ( i=wi-widthl; i<=wi+widthr; ++i ) {
		ii = f_domain(i, nx);
		jj = f_domain(j, ny);
		kk = f_domain(k, nz);
		pos[ii][jj][kk] = colour; 
		}
	      }
	    }
	  }
	  else {
	    for ( k=wk+vg_k; k<=wk; ++k ) {
	      for ( j=wj-widthl; j<=wj+widthr; ++j ) {
		for ( i=wi-widthl; i<=wi+widthr; ++i ) {
		  ii = f_domain(i, nx);
		  jj = f_domain(j, ny);
		  kk = f_domain(k, nz);
		  pos[ii][jj][kk] = colour; 
		}
	      }
	    }
	  }
	}
	else {
	  printf("Attention, wrong definition of the sqr_column %d\n",sqr_clm);
	  exit(1);
	}
      }
      else {
	printf("Attention, wrong definition of the sqr_column %d\n",sqr_clm);
	exit(1);
      }
      /* end of sqr_column */
      
    } else if ( strcmp(shape,"circ_column") == 0 ) { 
      fscanf(f0, "%f %f %f", &sx, &sy, &sz);
      fscanf(f0, "%d %d %d", &vg_i, &vg_j, &vg_k);
      fscanf(f0, "%f", &radious);
      fscanf(f0, "%f", &colour);
      circ_clm += 1;
      printf("circ_column %d: %f %f %f %d %d %d %f %f\n", circ_clm, sx, sy, sz, vg_i, vg_j, vg_k, radious, colour); 
      
      width = radious + 2;
      wi = sx;
      wj = sy;
      wk = sz;
      
      if ( vg_i!=0 && vg_j==0 && vg_k==0 ){
	if ( vg_i>0 ) {
	  for ( k=wk-width; k<=wk+width; ++k ) {
	    for ( j=wj-width; j<=wj+width; ++j ) {
	      for ( i=wi; i<=wi+vg_i; ++i ) {
		ii = f_domain(i, nx);
		jj = f_domain(j, ny);
		kk = f_domain(k, nz);
		if ( pow(radious,0.2e1) >= f_circ(sy, sz, jj, kk) ) {
		  pos[ii][jj][kk] = colour; 
		}
	      }
	    }
	  }
	}
	else {
	  for ( k=wk-width; k<=wk+width; ++k ) {
	    for ( j=wj-width; j<=wj+width; ++j ) {
	      for ( i=wi+vg_i; i<=wi; ++i ) {
		ii = f_domain(i, nx);
		jj = f_domain(j, ny);
		kk = f_domain(k, nz);
		if ( pow(radious,0.2e1) >= f_circ(sy, sz, jj, kk) ) {
		  pos[ii][jj][kk] = colour; 
		}
	      }
	    }
	  }
	}
      }
      else if ( vg_i==0 && vg_j!=0 && vg_k==0 ){
	if ( vg_j>0 ) {
	  for ( k=wk-width; k<=wk+width; ++k ) {
	    for ( j=wj; j<=wj+vg_j; ++j ) {
	      for ( i=wi-width; i<=wi+width; ++i ) {
		ii = f_domain(i, nx);
		jj = f_domain(j, ny);
		kk = f_domain(k, nz);
		if ( pow(radious,0.2e1) >= f_circ(sx, sz, ii, kk) ) {
		  pos[ii][jj][kk] = colour; 
		}
	      }
	    }
	  }
	}
	else {
	  for ( k=wk-width; k<=wk+width; ++k ) {
	    for ( j=wj+vg_j; j<=wj; ++j ) {
	      for ( i=wi-width; i<=wi+width; ++i ) {
		ii = f_domain(i, nx);
		jj = f_domain(j, ny);
		kk = f_domain(k, nz);
		if ( pow(radious,0.2e1) >= f_circ(sx, sz, ii, kk) ) {
		  pos[ii][jj][kk] = colour; 
		}
	      }
	    }
	  }
	}
      }
      else if ( vg_i==0 && vg_j==0 && vg_k!=0 ){
	if ( vg_k>0 ) {
	  for ( k=wk; k<=wk+vg_k; ++k ) {
	    for ( j=wj-width; j<=wj+width; ++j ) {
	      for ( i=wi-width; i<=wi+width; ++i ) {
		ii = f_domain(i, nx);
		jj = f_domain(j, ny);
		kk = f_domain(k, nz);
		if ( pow(radious,0.2e1) >= f_circ(sx, sy, ii, jj) ) {
		  pos[ii][jj][kk] = colour; 
		}
	      }
	    }
	  }
	}
	else {
	  for ( k=wk+vg_k; k<=wk; ++k ) {
	    for ( j=wj-width; j<=wj+width; ++j ) {
	      for ( i=wi-width; i<=wi+width; ++i ) {
		ii = f_domain(i, nx);
		jj = f_domain(j, ny);
		kk = f_domain(k, nz);
		if ( pow(radious,0.2e1) >= f_circ(sx, sy, ii, jj) ) {
		  pos[ii][jj][kk] = colour; 
		}
	      }
	    }
	  }
	}
      }
      else {
	printf("Attention, wrong definition of the sqr_column %d\n",sqr_clm);
	exit(1);
      }
    /* end of circ_column */
    }

  } 

  printf("all the geometries were readed and created\n");
  printf("number of wall: %d\n", wll);
  printf("number of cube_trj: %d\n", cb_trj);
  printf("number of rctgl_trj: %d\n", rctgl_trj);
  printf("number of sphr_trj: %d\n", sp_trj);
  printf("number of cube: %d\n", cb);
  printf("number of sphere: %d\n", sphr);
  printf("number of ellipse: %d\n", elps);
  printf("number of sqr_column: %d\n", sqr_clm);
  printf("number of circ_column: %d\n", circ_clm);

  
  fclose(f0);
  
  xdrstdio_create( &xdrs, f2, XDR_ENCODE ); 
  
  for ( k=0; k<nz; ++k ) {
    for (j=0; j<ny; ++j ) {
      for ( i=0; i<nx; ++i ) {
	/* Check output of xdr file*/
	xdr_float(&xdrs,&pos[i][j][k]); 
	fprintf(f1," %d %d %d %f \n", i,j,k,pos[i][j][k]);
      }     
    }
  }
  xdr_destroy(&xdrs); 
  
  fclose(f1);
  fclose(f2); 
  return 0;
}


/* function to keep the index  i, j, and k inside the computational domain */
int f_domain(int c1, int cn){
  int temp;
  temp=c1;
  if(c1 < 0){
   temp = 0;
  }
  if(c1 > cn-1){
    temp = cn-1;
  }
  return (temp);
}


/* function for 3D channel */
/*  x=x1+(x2-x1)/tmax*t 
    y=y1+(y2-y1)/tmax*t  
    z=z1+(z2-z1)/tmax*t */ 
int f_trj(int c1, int c2, int pmax, int p) {
  if( pmax == 0 ) {
    printf("Attention, wrong f_channel, there is no channel defined\n");
    exit(1);
  }
  else {
    return ( (float)c1 + (float)(c2-c1)*(float)p / (float)pmax );
  }
}


/* function for 3D sphere */
/* (sx-i)**2 + (sy-j)**2 + (sz-k)**2  */
float f_sphere(float cx, float cy, float cz, int i, int j, int k) {
  return (pow((float)cx-(float)i,0.2e1) + pow((float)cy-(float)j,0.2e1) + pow((float)cz-(float)k,0.2e1));
}
 
/* function for 3D circ */
/* (sx-i)**2 + (sy-j)**2  */
float f_circ(float cx, float cy, int i, int j) {
  return (pow((float)cx-(float)i,0.2e1) + pow((float)cy-(float)j,0.2e1));
}

/* function for 3D ellipse */
/*  maple routine */
/* 
with(CodeGeneration);

# the point distance (x1,y1,z1) and (x2,y2,z2) is constant 

dist:=((cx1-i)^2+(cy1-j)^2+(cz1-k)^2)^(1/2)+((cx2-i)^2+(cy2-j)^2+(cz2-k)^2)^(1/2);

*/
float f_ellipse(float cx1, float cy1, float cz1, float cx2, float cy2, float cz2, int i, int j, int k, float radious) {
  float temp1,temp2;
  if ( cx1 == cx2 && cy1 == cy2 && cz1 == cz2 ) {
    printf("Attention, wrong f_ellipse, (cx1,cy1,cz1)=(cx2,cy2,cz2)\n");
    exit(1);
  }
  else if (pow((float)cx1 -(float)cx2,0.2e1)+pow((float)cy1 -(float)cy2,0.2e1)+pow((float)cz1 -(float)cz2,0.2e1) > pow(radious,0.2e1) ) {
    printf("Attention, wrong f_ellipse, |(cx1,cy1,cz1)-(cx2,cy2,cz2)|< radious\n");
    exit(1);
  }
  else {
    temp1 = pow((float)cx1 - (float)i, 0.2e1) + pow((float)cy1 - (float)j, 0.2e1) + pow((float)cz1 - (float)k, 0.2e1);
    temp2 = pow((float)cx2 - (float)i, 0.2e1) + pow((float)cy2 - (float)j, 0.2e1) + pow((float)cz2 - (float)k, 0.2e1);
    return (sqrt(temp1)+sqrt(temp2));
  }
}




