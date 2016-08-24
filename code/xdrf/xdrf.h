/*_________________________________________________________________
 |
 | xdrf.h - include file for C routines that want to use the 
 |	    functions below.
*/

int xdropen(XDR *xdrs, const char *filename, const char *type);
int xdrclose(XDR *xdrs) ;
int xdr3dfcoord(XDR *xdrs, float *fp, int *size, float *precision) ;

