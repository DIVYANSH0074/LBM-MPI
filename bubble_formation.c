#include<stdio.h>
#include<mpi.h>
#include<stdlib.h>
#include<math.h>
#include<fstream>
#include<string>
#include<sstream>
#include<time.h>

int NX = 64;
int NY = 64;
int NZ = 64;
int ndir = 19;

double rd = 10.0;

double w0 = 1.0/3.0;
double ws = 1.0/18.0;
double wd = 1.0/36.0;

double G = -5.5;
double rhol = 2.0;
double rhov = 0.15;
double omega = 1.0;
double rho0 = 1.0;  

double tau = 1.0;
int tend = 1;

inline size_t scalar_index(int x,int y,int z)
{
 return NX*NY*z + NY*y + x;
}

inline size_t psi_index(int x,int y,int z)
{
 return NX*NY*(z+1) + NY*y + x;
}

inline size_t field_index(int x,int y,int z,int d)
{
 return ndir*(NX*NY*(z+1) + NX*y + x) + d;
}

void init_pdf(double *f,int rank_ny,int rank_ystart)
{
 for(int z=0,k=rank_ystart; z<rank_ny && k<NZ; z++,k++)
 {
  for(int y=0; y<NY; y++)
  {
   for(int x=0; x<NX;x++)
   {
   	double xr = x-NX/2;
   	double yr = y-NY/2;
    double zr = k-NZ/2;

    if(xr*xr + yr*yr + zr*zr <= rd*rd)
    {
     f[field_index(x,y,z,0)] = w0*rhol;

     f[field_index(x,y,z,1)] = ws*rhol;
     f[field_index(x,y,z,2)] = ws*rhol;
     f[field_index(x,y,z,3)] = ws*rhol;
     f[field_index(x,y,z,4)] = ws*rhol;
     f[field_index(x,y,z,5)] = ws*rhol;
     f[field_index(x,y,z,6)] = ws*rhol;

     f[field_index(x,y,z,7)] = wd*rhol;
     f[field_index(x,y,z,8)] = wd*rhol;
     f[field_index(x,y,z,9)] = wd*rhol;
     f[field_index(x,y,z,10)] = wd*rhol;
     f[field_index(x,y,z,11)] = wd*rhol;
     f[field_index(x,y,z,12)] = wd*rhol;
     f[field_index(x,y,z,13)] = wd*rhol;
     f[field_index(x,y,z,14)] = wd*rhol;
     f[field_index(x,y,z,15)] = wd*rhol;
     f[field_index(x,y,z,16)] = wd*rhol;
     f[field_index(x,y,z,17)] = wd*rhol;
     f[field_index(x,y,z,18)] = wd*rhol;
    }

    else
    {
     f[field_index(x,y,z,0)] = w0*rhov;

     f[field_index(x,y,z,1)] = ws*rhov;
     f[field_index(x,y,z,2)] = ws*rhov;
     f[field_index(x,y,z,3)] = ws*rhov;
     f[field_index(x,y,z,4)] = ws*rhov;
     f[field_index(x,y,z,5)] = ws*rhov;
     f[field_index(x,y,z,6)] = ws*rhov;

     f[field_index(x,y,z,7)] = wd*rhov;
     f[field_index(x,y,z,8)] = wd*rhov;
     f[field_index(x,y,z,9)] = wd*rhov;
     f[field_index(x,y,z,10)] = wd*rhov;
     f[field_index(x,y,z,11)] = wd*rhov;
     f[field_index(x,y,z,12)] = wd*rhov;
     f[field_index(x,y,z,13)] = wd*rhov;
     f[field_index(x,y,z,14)] = wd*rhov;
     f[field_index(x,y,z,15)] = wd*rhov;
     f[field_index(x,y,z,16)] = wd*rhov;
     f[field_index(x,y,z,17)] = wd*rhov;
     f[field_index(x,y,z,18)] = wd*rhov;
    }
   }
  }
 }
}

void stream_save(double *f1,double *rho,double *ux,double *uy,double *uz,double *psi,int rank_ny)
{
 for(int z=0; z<rank_ny; z++)
 {
  for(int y=0; y<NY; y++)
  {
   for(int x=0; x<NX; x++)
   {
    int xm = (NX+x-1)%NX;
    int xp = (x+1)%NX;

    int ym = (NY+y-1)%NY;
    int yp = (y+1)%NY;

    int zm = (z-1)%NZ;
    int zp = (z+1)%NZ;

    double ft0 = f1[field_index(x,y,z,0)];

    double ft1 = f1[field_index(xm,y,z,1)];
    double ft2 = f1[field_index(x,ym,z,2)];
    double ft3 = f1[field_index(x,y,zm,3)];
    double ft4 = f1[field_index(xp,y,z,4)];
    double ft5 = f1[field_index(x,yp,z,5)];
    double ft6 = f1[field_index(x,y,zp,6)];

    double ft7 = f1[field_index(xm,ym,z,7)];
    double ft8 = f1[field_index(xp,ym,z,8)];
    double ft9 = f1[field_index(xp,yp,z,9)];
    double ft10 = f1[field_index(xm,yp,z,10)];

    double ft11 = f1[field_index(x,ym,zm,11)];
    double ft12 = f1[field_index(x,yp,zm,12)];
    double ft13 = f1[field_index(x,yp,zp,13)];
    double ft14 = f1[field_index(x,ym,zp,14)];

    double ft15 = f1[field_index(xm,y,zm,15)];
    double ft16 = f1[field_index(xp,y,zm,16)];
    double ft17 = f1[field_index(xp,y,zp,17)];
    double ft18 = f1[field_index(xm,y,zp,18)];

    double r = ft0+ft1+ft2+ft3+ft4+ft5+ft6+ft7+ft8+ft9+ft10+ft11+ft12+ft13+ft14+ft15+ft16+ft17+ft18;
    double rinv = 1.0/r;
    double u = rinv*(ft1+ft7+ft10+ft15+ft18 - (ft4+ft8+ft9+ft16+ft17));
    double v = rinv*(ft2+ft7+ft8+ft11+ft14 - (ft5+ft9+ft10+ft12+ft13));
    double w = rinv*(ft3+ft11+ft12+ft15+ft16 - (ft6+ft13+ft14+ft17+ft18));

    rho[scalar_index(x,y,z)] = r;
    ux[scalar_index(x,y,z)] = u;
    uy[scalar_index(x,y,z)] = v;
    uz[scalar_index(x,y,z)] = w;

    double p = rho0*(1.0-exp(-r/rho0));
    
    psi[psi_index(x,y,z)] = p;
   }
  }
 }
}

void collide(double *f1,double *f2,double *rho,double *ux,double *uy,double *uz,double *psi,int rank_ny)
{
 for(int z=0; z<rank_ny; z++)
 {
  for(int y=0; y<NY; y++)
  {
   for(int x=0; x<NX; x++)
   {
    int xm = (NX+x-1)%NX;
    int xp = (x+1)%NX;

    int ym = (NY+y-1)%NY;
    int yp = (y+1)%NY;

    int zm = (z-1)%NZ;
    int zp = (z+1)%NZ;

    double ft0 = f1[field_index(x,y,z,0)];

    double ft1 = f1[field_index(xm,y,z,1)];
    double ft2 = f1[field_index(x,ym,z,2)];
    double ft3 = f1[field_index(x,y,zm,3)];
    double ft4 = f1[field_index(xp,y,z,4)];
    double ft5 = f1[field_index(x,yp,z,5)];
    double ft6 = f1[field_index(x,y,zp,6)];

    double ft7 = f1[field_index(xm,ym,z,7)];
    double ft8 = f1[field_index(xp,ym,z,8)];
    double ft9 = f1[field_index(xp,yp,z,9)];
    double ft10 = f1[field_index(xm,yp,z,10)];

    double ft11 = f1[field_index(x,ym,zm,11)];
    double ft12 = f1[field_index(x,yp,zm,12)];
    double ft13 = f1[field_index(x,yp,zp,13)];
    double ft14 = f1[field_index(x,ym,zp,14)];

    double ft15 = f1[field_index(xm,y,zm,15)];
    double ft16 = f1[field_index(xp,y,zm,16)];
    double ft17 = f1[field_index(xp,y,zp,17)];
    double ft18 = f1[field_index(xm,y,zp,18)];


    double ps0 = psi[psi_index(x,y,z)];

    double ps1 = psi[psi_index(xp,y,z)];
    double ps2 = psi[psi_index(x,yp,z)];
    double ps3 = psi[psi_index(x,y,zp)];
    double ps4 = psi[psi_index(xm,y,z)];
    double ps5 = psi[psi_index(x,ym,z)];
    double ps6 = psi[psi_index(x,y,zm)];

    double ps7 = psi[psi_index(xp,yp,z)];
    double ps8 = psi[psi_index(xm,yp,z)];
    double ps9 = psi[psi_index(xm,ym,z)];
    double ps10 = psi[psi_index(xp,ym,z)];

    double ps11 = psi[psi_index(x,yp,zp)];
    double ps12 = psi[psi_index(x,ym,zp)];
    double ps13 = psi[psi_index(x,ym,zm)];
    double ps14 = psi[psi_index(x,yp,zm)];

    double ps15 = psi[psi_index(xp,y,zp)];
    double ps16 = psi[psi_index(xm,y,zp)];
    double ps17 = psi[psi_index(xm,y,zm)];
    double ps18 = psi[psi_index(xp,y,zm)];

    double r = rho[scalar_index(x,y,z)];
    double u = ux[scalar_index(x,y,z)];
    double v = uy[scalar_index(x,y,z)];
    double w = uz[scalar_index(x,y,z)];

    double fx = ws*(ps1-ps4) + wd*(ps7+ps10+ps15+ps18 - (ps8+ps9+ps16+ps17));
    double fy = ws*(ps2-ps5) + wd*(ps7+ps8+ps11+ps14 - (ps9+ps10+ps12+ps13));
    double fz = ws*(ps3-ps6) + wd*(ps11+ps12+ps15+ps16 - (ps13+ps14+ps17+ps18));

    fx = -G*fx*ps0;
    fy = -G*fy*ps0;
    fz = -G*fz*ps0;

    double ueq = u + fx*tau/r;
    double veq = v + fy*tau/r;
    double weq = w + fz*tau/r;

    double tw0r = omega*w0*r;
    double twsr = omega*ws*r;
    double twdr = omega*wd*r;

    double omusq = 1.0-1.5*(ueq*ueq+veq*veq+weq*weq);
   
    double tux = 3.0*ueq;
    double tuy = 3.0*veq;
    double tuz = 3.0*weq;

    double om = 1.0 - omega;

    f2[field_index(x,y,z,0)] = om*ft0 + tw0r*omusq;

    double cidot3u = tux;
    f2[field_index(x,y,z,1)] = om*ft1 + twsr*(omusq + cidot3u*(1.0 + 0.5*cidot3u));
    f2[field_index(x,y,z,4)] = om*ft4 + twsr*(omusq - cidot3u*(1.0 - 0.5*cidot3u));

    cidot3u = tuy;
    f2[field_index(x,y,z,2)] = om*ft2 + twsr*(omusq + cidot3u*(1.0 + 0.5*cidot3u));
    f2[field_index(x,y,z,5)] = om*ft5 + twsr*(omusq - cidot3u*(1.0 - 0.5*cidot3u));

    cidot3u = tuz;
    f2[field_index(x,y,z,3)] = om*ft3 + twsr*(omusq + cidot3u*(1.0 + 0.5*cidot3u));
    f2[field_index(x,y,z,6)] = om*ft6 + twsr*(omusq - cidot3u*(1.0 - 0.5*cidot3u));

    cidot3u = tux+tuy;
    f2[field_index(x,y,z,7)] = om*ft7 + twdr*(omusq + cidot3u*(1.0 + 0.5*cidot3u));
    f2[field_index(x,y,z,9)] = om*ft9 + twdr*(omusq - cidot3u*(1.0 - 0.5*cidot3u));

    cidot3u = -tux+tuy;
    f2[field_index(x,y,z,8)] = om*ft8 + twdr*(omusq + cidot3u*(1.0 + 0.5*cidot3u));
    f2[field_index(x,y,z,10)] = om*ft10 + twdr*(omusq - cidot3u*(1.0 - 0.5*cidot3u));

    cidot3u = tux+tuz;
    f2[field_index(x,y,z,15)] = om*ft15 + twdr*(omusq + cidot3u*(1.0 + 0.5*cidot3u));
    f2[field_index(x,y,z,17)] = om*ft17 + twdr*(omusq - cidot3u*(1.0 - 0.5*cidot3u));

    cidot3u = -tux+tuz;
    f2[field_index(x,y,z,16)] = om*ft16 + twdr*(omusq + cidot3u*(1.0 + 0.5*cidot3u));
    f2[field_index(x,y,z,18)] = om*ft18 + twdr*(omusq - cidot3u*(1.0 - 0.5*cidot3u));

    cidot3u = tuy+tuz;
    f2[field_index(x,y,z,11)] = om*ft11 + twdr*(omusq + cidot3u*(1.0 + 0.5*cidot3u));
    f2[field_index(x,y,z,13)] = om*ft13 + twdr*(omusq - cidot3u*(1.0 - 0.5*cidot3u));

    cidot3u = -tuy+tuz;
    f2[field_index(x,y,z,12)] = om*ft12 + twdr*(omusq + cidot3u*(1.0 + 0.5*cidot3u));
    f2[field_index(x,y,z,14)] = om*ft14 + twdr*(omusq - cidot3u*(1.0 - 0.5*cidot3u));
   }
  }
 }
}

double* gather(double *rho,int rank_ny,int rank_ystart,int size,int rank)
{
 int tag;
 MPI_Status status;
 double *r = NULL;
 int *rank_ny0 = NULL,*rank_ystart0 = NULL;

 if(rank==0)
 {
  r = (double*) malloc(NX*NY*NZ*sizeof(double));
  rank_ny0 = (int*) malloc(size*sizeof(int));
  rank_ystart0 = (int*) malloc(size*sizeof(int));
 }

 if(rank>0)
 {
  tag = 1;
  MPI_Send(&rank_ny,1,MPI_INT,0,tag,MPI_COMM_WORLD);
  tag = 2;
  MPI_Send(&rank_ystart,1,MPI_INT,0,tag,MPI_COMM_WORLD);
  tag = 3;
  MPI_Send(&rho[scalar_index(0,0,0)],NX*NY*rank_ny,MPI_DOUBLE,0,tag,MPI_COMM_WORLD);
 }
 
 if(rank==0)
 {
  tag = 1;
  for(int i=1; i<size; i++)
  {
   MPI_Recv(&rank_ny0[i],1,MPI_INT,i,tag,MPI_COMM_WORLD,&status);
  }
  tag = 2;
  for(int i=1; i<size; i++)
  {
   MPI_Recv(&rank_ystart0[i],1,MPI_INT,i,tag,MPI_COMM_WORLD,&status);
  }
 }

 if(rank==0)
 {
  tag = 3;
  for(int i=1; i<size; i++)
  {
   MPI_Recv(&r[NX*NY*rank_ystart0[i]],NX*NY*rank_ny0[i],MPI_DOUBLE,i,tag,MPI_COMM_WORLD,&status);
  }
  
  for(int z=0; z<rank_ny; z++)
  {
   for(int y=0; y<NY; y++)
   {
    for(int x=0; x<NX; x++)
    {
     r[scalar_index(x,y,z)] = rho[scalar_index(x,y,z)];
    }
   }
  }

  return r;
 }
 
 free(r);
 free(rank_ystart0);
 free(rank_ny0);
}

void SaveVTKImageData_ascii(double *rho1,std::string fileName,std::string dataName,float *origin,int *spacing,int *dims)
{
  FILE* pFile;
  pFile = fopen(fileName.c_str(), "w");
  fprintf(pFile, "# vtk DataFile Version 3.0 \n");
  fprintf(pFile, "Slow ASCII version...\n");
  fprintf(pFile, "ASCII \n");
  fprintf(pFile, "\n");
  fprintf(pFile, "DATASET STRUCTURED_POINTS \n");
  fprintf(pFile, "DIMENSIONS %d  %d  %d \n", dims[0], dims[1], dims[2]);
  fprintf(pFile, "\n");
  fprintf(pFile, "ORIGIN  %4.3f   %4.3f  %4.3f \n", origin[0], origin[1], origin[2]);
  fprintf(pFile, "SPACING  %d  %d  %d \n", spacing[0], spacing[1], spacing[2]);
  fprintf(pFile, "\n");
  fprintf(pFile, "POINT_DATA %d  \n", dims[0] * dims[1] * dims[2]);
  fprintf(pFile, "SCALARS \t %s  \t float \n", dataName.c_str());
  fprintf(pFile, "LOOKUP_TABLE default \n \n");

  for (int z = 0; z < NZ; z++) {
    for (int y = 0; y < NY; y++) {
      for (int x = 0; x < NX; x++) {

        fprintf(pFile, " %g  ", rho1[scalar_index(x, y, z)]);
      }
    }
    fprintf(pFile, "\n");
  }

  fclose(pFile);
}

int main(int argc,char* argv[])
{
 int rank,size;
 MPI_Init(&argc,&argv);
 MPI_Comm_size(MPI_COMM_WORLD,&size);
 MPI_Comm_rank(MPI_COMM_WORLD,&rank);

 int rank_ny,rank_ystart;

 if(rank < NZ % size)
 {
  rank_ny = NZ/size + 1;
  rank_ystart = rank * rank_ny;
 }
 
 else
 {
  rank_ny = NZ/size;
  rank_ystart = NZ - (size-rank)*rank_ny;
 }

 double *f1 = (double*) malloc(NX*NY*(rank_ny+2)*ndir*sizeof(double));
 double *f2 = (double*) malloc(NX*NY*(rank_ny+2)*ndir*sizeof(double));
 double *rho = (double*) malloc(NX*NY*rank_ny*sizeof(double));
 double *ux = (double*) malloc(NX*NY*rank_ny*sizeof(double));
 double *uy = (double*) malloc(NX*NY*rank_ny*sizeof(double));
 double *uz = (double*) malloc(NX*NY*rank_ny*sizeof(double));
 double *psi = (double*) malloc(NX*NY*(rank_ny+2)*sizeof(double));
 double *rhof = (double*) malloc(NX*NY*NZ*sizeof(double));

 int rankp1 = (rank+1) % size;
 int rankm1 = (size+rank-1) % size;

 init_pdf(f1,rank_ny,rank_ystart);

 size_t transfer_doubles = ndir*NX*NY;
 size_t transfer_psi = NX*NY;

 clock_t start, end;
 double time;

 start = clock();

 for(int t=0; t<tend; t++)
 {

  MPI_Sendrecv(&f1[field_index(0,0,rank_ny-1,0)],transfer_doubles,MPI_DOUBLE,rankp1,rank,&f1[field_index(0,0,-1,0)],
  transfer_doubles,MPI_DOUBLE,rankm1,rankm1,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
  
  MPI_Sendrecv(&f1[field_index(0,0,0,0)],transfer_doubles,MPI_DOUBLE,rankm1,rank,&f1[field_index(0,0,rank_ny,0)],
  transfer_doubles,MPI_DOUBLE,rankp1,rankp1,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
 
  stream_save(f1,rho,ux,uy,uz,psi,rank_ny);

  MPI_Sendrecv(&psi[psi_index(0,0,rank_ny-1)],transfer_psi,MPI_DOUBLE,rankp1,rank,&psi[psi_index(0,0,-1)],
  transfer_psi,MPI_DOUBLE,rankm1,rankm1,MPI_COMM_WORLD,MPI_STATUS_IGNORE);

  MPI_Sendrecv(&psi[psi_index(0,0,0)],transfer_psi,MPI_DOUBLE,rankm1,rank,&psi[psi_index(0,0,rank_ny)],
  transfer_psi,MPI_DOUBLE,rankp1,rankp1,MPI_COMM_WORLD,MPI_STATUS_IGNORE);

  collide(f1,f2,rho,ux,uy,uz,psi,rank_ny);

  double* temp = f1;
  f1 = f2;
  f2 = temp;
 }

 end = clock();

 time = ((double) (end-start)) / CLOCKS_PER_SEC;
 printf("%f\n",time);

 rhof = gather(rho,rank_ny,rank_ystart,size,rank); 

 if(rank==0)
 {

  std::string densityFileStub("Density");
  std::string fileSuffix(".vtk");
  std::stringstream ts_ind;
  std::string ts_ind_str;
  int l_conv_fact = 1;
  int vtk_ts = 0;
  std::string fileName;
  std::string dataName("density");
  int dims[3];
  dims[0] = NX; dims[1] = NY; dims[2] = NZ;
  float origin[3];
  origin[0] = 0.; origin[1] = 0.; origin[2] = 0.;
  int spacing[3];
  spacing[0] = l_conv_fact; spacing[1] = l_conv_fact; spacing[2] = l_conv_fact;

 //Writing VTK File
  ts_ind << vtk_ts;

  fileName = densityFileStub + ts_ind.str() + fileSuffix;
  ts_ind.str("");
  SaveVTKImageData_ascii(rhof, fileName, dataName, origin, spacing, dims);

 } 

 free(f1);
 free(f2);
 free(rho);
 free(ux);
 free(uy);
 free(uz);
 free(psi);
 free(rhof);

 MPI_Finalize();
 return 0;
}