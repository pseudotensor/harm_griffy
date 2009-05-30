/* 
	produces binary files.
*/

#ifdef DEBUG
#define ZLOOP_DUMP ZLOOPG
#else
#define ZLOOP_DUMP ZLOOP
#endif

#include "decs.h"
void bdump_field(const char *binnam, double field[INA][JNA])
{
  int i,j ;
  static size_t size=sizeof(double);
  FILE *binary_file;

  binary_file = fopen(binnam,"wb") ;
  if(binary_file==NULL) {
    fprintf(stderr,"error opening image file\n") ;
    exit(END_BINARY_OUTPUT_OPENING_FILE) ;
  }
  ZLOOP_DUMP { fwrite(&(field[i][j]),size,1,binary_file) ; }
  fclose(binary_file) ;
}

void bdump_primitive(const char *binnam,const int field)
{
  int i,j ;
  static size_t size=sizeof(p[IS][IS][IS]);
  FILE *binary_file;

  binary_file = fopen(binnam,"wb") ;
  if(binary_file==NULL) {
    fprintf(stderr,"error opening image file\n") ;
    exit(END_BINARY_OUTPUT_OPENING_FILE) ;
  }
  ZLOOP_DUMP { fwrite(&(p[i][j][field]),size,1,binary_file) ; }
  fclose(binary_file) ;
}

void bdump_conn(const char *binnam,
		const int idim,const int jdim,const int kdim)
{
  int i,j ;
  static size_t size=sizeof(p[IS][IS][IS]);
  FILE *binary_file;

  binary_file = fopen(binnam,"wb") ;
  if(binary_file==NULL) {
    fprintf(stderr,"error opening image file\n") ;
    exit(END_BINARY_OUTPUT_OPENING_FILE) ;
  }
  ZLOOP_DUMP { fwrite(&(conn[i][j][idim][jdim][kdim]),size,1,binary_file) ; }
  fclose(binary_file) ;
}

void bdump_gcon(const char *binnam,const int loc,const int jdim,const int kdim)
{
  int i,j ;
  static size_t size=sizeof(p[IS][IS][IS]);
  FILE *binary_file;

  binary_file = fopen(binnam,"wb") ;
  if(binary_file==NULL) {
    fprintf(stderr,"error opening image file\n") ;
    exit(END_BINARY_OUTPUT_OPENING_FILE) ;
  }
  ZLOOP_DUMP { fwrite(&(gcon[i][j][loc][jdim][kdim]),size,1,binary_file) ; }
  fclose(binary_file) ;
}

void bdump_gcov(const char *binnam,const int loc,const int jdim,const int kdim)
{
  int i,j ;
  static size_t size=sizeof(p[IS][IS][IS]);
  FILE *binary_file;

  binary_file = fopen(binnam,"wb") ;
  if(binary_file==NULL) {
    fprintf(stderr,"error opening image file\n") ;
    exit(END_BINARY_OUTPUT_OPENING_FILE) ;
  }
  ZLOOP_DUMP { fwrite(&(gcov[i][j][loc][jdim][kdim]),size,1,binary_file) ; }
  fclose(binary_file) ;
}

void bdump_gdet(const char *binnam,const int loc)
{
  int i,j ;
  static size_t size=sizeof(p[IS][IS][IS]);
  FILE *binary_file;

  binary_file = fopen(binnam,"wb") ;
  if(binary_file==NULL) {
    fprintf(stderr,"error opening image file\n") ;
    exit(END_BINARY_OUTPUT_OPENING_FILE) ;
  }
  ZLOOP_DUMP { fwrite(&(gdet[i][j][loc]),size,1,binary_file) ; }
  fclose(binary_file) ;
}

void bdump_array(const char *binnam,
		 const int n, double a[n],
		 const int begin, const int end){
  int index ;
  static size_t size=sizeof(double) ;
  FILE *binary_file;

  binary_file = fopen(binnam,"wb") ;
  if(binary_file==NULL) {
    fprintf(stderr,"error opening image file\n") ;
    exit(END_BINARY_OUTPUT_OPENING_FILE) ;
  }
  for(index=begin;index<=end;index++) {
    fwrite(&(a[index]),size,1,binary_file) ;
  }
  fclose(binary_file) ;
}

void bdump_fluxes(const int label)
{
#define FILENAMESIZE 100
    char binnam[FILENAMESIZE];
#undef FILENAMESIZE
    int i,j,k;
    FILE *binary_file;
    static size_t size=sizeof(double);
    PLOOP {
      sprintf(binnam,"b_F1_%01d__%02d_%04d",k,label,nstep);
      binary_file = fopen(binnam,"wb") ;
      ZLOOP_DUMP { fwrite(&(F1[i][j][k]),size,1,binary_file) ; }
      fclose(binary_file) ;
    }
    PLOOP {
      sprintf(binnam,"b_F2_%01d__%02d_%04d",k,label,nstep);
      binary_file = fopen(binnam,"wb") ;
      ZLOOP_DUMP { fwrite(&(F2[i][j][k]),size,1,binary_file) ; }
      fclose(binary_file) ;
    }
}

void bdump_source(const double dU_debug[INA][JNA][NPR] )
{
#define FILENAMESIZE 100
    char binnam[FILENAMESIZE];
#undef FILENAMESIZE
    int i,j,k;
    FILE *binary_file;
    static size_t size=sizeof(double);
    PLOOP {
      sprintf(binnam,"b_dU_%01d___%04d",k,nstep);
      binary_file = fopen(binnam,"wb") ;
      ZLOOP_DUMP { fwrite(&(dU_debug[i][j][k]),size,1,binary_file) ; }
      fclose(binary_file) ;
    }
}
