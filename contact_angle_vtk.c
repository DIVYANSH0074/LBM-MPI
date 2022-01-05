#include<stdio.h>
#include<mpi.h>
#include<stdlib.h>
#include<math.h>
#include<fstream>
#include<string>
#include<sstream>

const int NX = 50;
const int NY = 50;
const int NZ = 50;
const int ndir = 19;

const int solid_size = 10;

const double w0 = 1.0/3.0;
const double ws = 1.0/18.0;
const double wd = 1.0/36.0;

const double wi[] = {w0, ws, ws, ws, ws, ws, ws, wd, wd, wd, wd, wd, wd, wd, wd, wd, wd, wd, wd};

const int dirx[] = {0, 1, -1, 0, 0, 0, 0, 1, 1, 1, 1, -1, -1, -1, -1, 0, 0, 0, 0};
const int diry[] = {0, 0, 0, 1, -1, 0, 0, 1, -1, 0, 0, 1, -1, 0, 0, 1, 1, -1, -1};
const int dirz[] = {0, 0, 0, 0, 0, 1, -1, 0, 0, 1, -1, 0, 0, 1, -1, 1, -1, 1, -1};

const int oopf[] = {0,2,1,4,3,6,5,12,11,14,13,8,7,10,9,18,17,16,15};

const double G = 0.85;
const double Gads1 = -2;
const double Gads2 = 2;
const double rhol = 0.06;
const double rhov = 2.0;
const double omega = 1.0;

const double tau = 1.0;
const int tend = 200;

inline size_t scalar_index(int x,int y,int z)
{
 return NX*NY*z + NX*y + x;
}

inline size_t field_index(int x,int y,int z,int d)
{
 return ndir*(NX*NY*(z+1) + NX*y + x) + d;
} 

inline size_t sn_index(int x,int y,int z)
{
 return NX*NY*(z+1) + NX*y + x; 
}

void init_solid(int *solid,int rank_ny,int rank_ystart,int rank)
{
  for(int z=0,k=rank_ystart; z<rank_ny && k<NZ; z++,k++)
  {
   for(int y=0; y<NY; y++)
   {
    for(int x=0; x<NX; x++)
    {
     if(k<10)
     {
      solid[sn_index(x,y,z)] = 1;
     }

     else
     {
      solid[sn_index(x,y,z)] = 0;
     }
    }
   }
  }
}

void init_rho(double *rho1,double *rho2,int rank_ny,int rank_ystart)
{
 for(int z=0,k=rank_ystart; z<rank_ny && k<NZ; z++,k++)
 {
  for(int y=0; y<NY; y++)
  {
   for(int x=0; x<NX; x++)
   {
    if(k<20 && k>=solid_size)
    {
     if((x<NX/2+10 && x>=NX/2-11) && (y<NY/2+10 && y>=NY/2-11))
     {
      rho1[sn_index(x,y,z)] = rhov;
      rho2[sn_index(x,y,z)] = rhol;
     }

     else
     {
      rho1[sn_index(x,y,z)] = rhol;
      rho2[sn_index(x,y,z)] = rhov; 
     }
    }

    else
    {
     rho1[sn_index(x,y,z)] = rhol;
     rho2[sn_index(x,y,z)] = rhov;
    }
   }
  }
 }
}

void init_pdf(double *f1,double *f2,double *rho1,double *rho2,int rank_ny)
{
 for(int z=0; z<rank_ny; z++)
 {
  for(int y=0; y<NY; y++)
  {
   for(int x=0; x<NX; x++)
   {
    for(int d=0; d<ndir; d++)
    {
     f1[field_index(x,y,z,d)] = wi[d]*rho1[sn_index(x,y,z)];
     f2[field_index(x,y,z,d)] = wi[d]*rho2[sn_index(x,y,z)];
    }
   }
  }
 }
}

void stream(double *f1,double *fc1,double *f2,double *fc2,int rank_ny)
{
 for(int z=0; z<rank_ny; z++)
 {
  for(int y=0; y<NY; y++)
  {
   for(int x=0; x<NX; x++)
   {
    for(int d=0; d<ndir; d++)
    {
     int xmd = (NX+x-dirx[d])%NX;
     int ymd = (NY+y-diry[d])%NY;
     int zmd = (z-dirz[d]);

     fc1[field_index(x,y,z,d)] = f1[field_index(xmd,ymd,zmd,d)];

     fc2[field_index(x,y,z,d)] = f2[field_index(xmd,ymd,zmd,d)];
    }
   }
  }
 }
}

void compute_rho_u(double *fc1,double *fc2,double *rho1,double *rho2,double *ux,double *uy,double *uz,int rank_ny)
{
 for(int z=0; z<rank_ny; z++)
 {
  for(int y=0; y<NY; y++)
  {
   for(int x=0; x<NX; x++)
   {
    double r1 = 0.0;
    double r2 = 0.0;
    double u = 0.0;
    double v = 0.0;
    double w = 0.0;

    for(int d=0; d<ndir; d++)
    {
     r1 += fc1[field_index(x,y,z,d)];
     r2 += fc2[field_index(x,y,z,d)];

     u += dirx[d]*(fc1[field_index(x,y,z,d)]+fc2[field_index(x,y,z,d)]);
     v += diry[d]*(fc1[field_index(x,y,z,d)]+fc2[field_index(x,y,z,d)]);
     w += dirz[d]*(fc1[field_index(x,y,z,d)]+fc2[field_index(x,y,z,d)]); 
    }

    rho1[sn_index(x,y,z)] = r1;
    rho2[sn_index(x,y,z)] = r2;

    ux[scalar_index(x,y,z)] = u/(r1+r2);
    uy[scalar_index(x,y,z)] = v/(r1+r2);
    uz[scalar_index(x,y,z)] = w/(r1+r2);
   }
  }
 }
}

void collide(double *fc1,double *fc2,double *rho1,double *rho2,double *ux,double *uy,double *uz,int *solid,int rank_ny,int rank_ystart)
{
 double om = 1.0 - omega;

 for(int z=0,k=rank_ystart; z<rank_ny && k<NZ; z++,k++)
 {
  for(int y=0; y<NY; y++)
  {
   for(int x=0; x<NX; x++)
   {
    double r1 = rho1[sn_index(x,y,z)];
    double r2 = rho2[sn_index(x,y,z)];

    double u = ux[scalar_index(x,y,z)];
    double v = uy[scalar_index(x,y,z)];
    double w = uz[scalar_index(x,y,z)];

    double fx1 = 0.0;
    double fx2 = 0.0;

    double fy1 = 0.0;
    double fy2 = 0.0;

    double fz1 = 0.0;
    double fz2 = 0.0;

    double fsx1 = 0.0;
    double fsx2 = 0.0;

    double fsy1 = 0.0;
    double fsy2 = 0.0;

    double fsz1 = 0.0;
    double fsz2 = 0.0;
    
    for(int d=0; d<ndir; d++)
    {
     int xm = (NX+x+dirx[d])%NX;
     int ym = (NY+y+diry[d])%NY;
     int zm = (z+dirz[d])%NZ;

     fx1 += wi[d]*dirx[d]*rho2[sn_index(xm,ym,zm)];
     fx2 += wi[d]*dirx[d]*rho1[sn_index(xm,ym,zm)];

     fy1 += wi[d]*diry[d]*rho2[sn_index(xm,ym,zm)];
     fy2 += wi[d]*diry[d]*rho1[sn_index(xm,ym,zm)];

     fz1 += wi[d]*dirz[d]*rho2[sn_index(xm,ym,zm)];
     fz2 += wi[d]*dirz[d]*rho1[sn_index(xm,ym,zm)];

     fsx1 += wi[d]*dirx[d]*solid[sn_index(xm,ym,zm)];
     fsx2 += wi[d]*dirx[d]*solid[sn_index(xm,ym,zm)];

     fsy1 += wi[d]*diry[d]*solid[sn_index(xm,ym,zm)];
     fsy2 += wi[d]*diry[d]*solid[sn_index(xm,ym,zm)];

     fsz1 += wi[d]*dirz[d]*solid[sn_index(xm,ym,zm)];
     fsz2 += wi[d]*dirz[d]*solid[sn_index(xm,ym,zm)];
    }

    fx1 = -G*fx1*rho1[sn_index(x,y,z)];
    fx2 = -G*fx2*rho2[sn_index(x,y,z)];

    fy1 = -G*fy1*rho1[sn_index(x,y,z)];
    fy2 = -G*fy2*rho2[sn_index(x,y,z)];

    fz1 = -G*fz1*rho1[sn_index(x,y,z)];
    fz2 = -G*fz2*rho2[sn_index(x,y,z)];

    fsx1 = -Gads1*fsx1*rho1[sn_index(x,y,z)];
    fsx2 = -Gads2*fsx2*rho2[sn_index(x,y,z)];

    fsy1 = -Gads1*fsy1*rho1[sn_index(x,y,z)];
    fsy2 = -Gads2*fsy2*rho2[sn_index(x,y,z)];

    fsz1 = -Gads1*fsz1*rho1[sn_index(x,y,z)];
    fsz2 = -Gads2*fsz2*rho2[sn_index(x,y,z)];

    double ueq1 = u + (fx1+fsx1)*tau/r1;
    double ueq2 = u + (fx2+fsx2)*tau/r2;

    double veq1 = v + (fy1+fsy1)*tau/r1;
    double veq2 = v + (fy2+fsy2)*tau/r2;

    double weq1 = w + (fz1+fsz1)*tau/r1;
    double weq2 = w + (fz2+fsz2)*tau/r2;

    if(k<solid_size)
    {
     double aa[19] , bb[19];

     for(int d=0; d<ndir; d++)
     {
      aa[d] = fc1[field_index(x,y,z,oopf[d])];
      bb[d] = fc2[field_index(x,y,z,oopf[d])];
     }

     for(int d=0; d<ndir; d++)
     {
      fc1[field_index(x,y,z,d)] = aa[d];
      fc2[field_index(x,y,z,d)] = bb[d];
     }
    }

    else
    {
     for(int d=0; d<ndir; d++)
     {
      double cidotu1 = dirx[d]*ueq1 + diry[d]*veq1 + dirz[d]*weq1;
      double cidotu2 = dirx[d]*ueq2 + diry[d]*veq2 + dirz[d]*weq2;

      double feq1 = wi[d]*r1*(1.0 + 3.0*cidotu1 + 4.5*cidotu1*cidotu1 -1.5*(ueq1*ueq1+veq1*veq1+weq1*weq1));
      double feq2 = wi[d]*r2*(1.0 + 3.0*cidotu2 + 4.5*cidotu2*cidotu2 -1.5*(ueq2*ueq2+veq2*veq2+weq2*weq2));

      fc1[field_index(x,y,z,d)] = om*fc1[field_index(x,y,z,d)] + omega*feq1;
      fc2[field_index(x,y,z,d)] = om*fc2[field_index(x,y,z,d)] + omega*feq2;
     }
    }
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
  MPI_Send(&rho[sn_index(0,0,0)],NX*NY*rank_ny,MPI_DOUBLE,0,tag,MPI_COMM_WORLD);
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
     r[scalar_index(x,y,z)] = rho[sn_index(x,y,z)];
    }
   }
  }
 }

 return r;
 
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

 double *fc1 = (double*) malloc(NX*NY*(rank_ny+2)*ndir*sizeof(double));
 double *fc2 = (double*) malloc(NX*NY*(rank_ny+2)*ndir*sizeof(double));

 int *solid = (int*) malloc(NX*NY*(rank_ny+2)*sizeof(int));

 double *rho1 = (double*) malloc(NX*NY*(rank_ny+2)*sizeof(double));
 double *rho2 = (double*) malloc(NX*NY*(rank_ny+2)*sizeof(double));
 double *rhof = (double*) malloc(NX*NY*NZ*sizeof(double));

 double *ux = (double*) malloc(NX*NY*rank_ny*sizeof(double));
 double *uy = (double*) malloc(NX*NY*rank_ny*sizeof(double));
 double *uz = (double*) malloc(NX*NY*rank_ny*sizeof(double));

 int rankp1 = (rank+1) % size;
 int rankm1 = (size+rank-1) % size;

 init_solid(solid,rank_ny,rank_ystart,rank);

 init_rho(rho1,rho2,rank_ny,rank_ystart);

 init_pdf(f1,f2,rho1,rho2,rank_ny);

 size_t transfer_doubles = ndir*NX*NY;
 size_t transfer_prop = NX*NY;

 MPI_Sendrecv(&solid[sn_index(0,0,rank_ny-1)],transfer_prop,MPI_INT,rankp1,rank,&solid[sn_index(0,0,-1)],
 transfer_prop,MPI_INT,rankm1,rankm1,MPI_COMM_WORLD,MPI_STATUS_IGNORE);

 MPI_Sendrecv(&solid[sn_index(0,0,0)],transfer_prop,MPI_INT,rankm1,rank,&solid[sn_index(0,0,rank_ny)],
 transfer_prop,MPI_INT,rankp1,rankp1,MPI_COMM_WORLD,MPI_STATUS_IGNORE);

 for(int t=0; t<tend; t++)
 {

  MPI_Sendrecv(&f1[field_index(0,0,rank_ny-1,0)],transfer_doubles,MPI_DOUBLE,rankp1,rank,&f1[field_index(0,0,-1,0)],
  transfer_doubles,MPI_DOUBLE,rankm1,rankm1,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
  
  MPI_Sendrecv(&f1[field_index(0,0,0,0)],transfer_doubles,MPI_DOUBLE,rankm1,rank,&f1[field_index(0,0,rank_ny,0)],
  transfer_doubles,MPI_DOUBLE,rankp1,rankp1,MPI_COMM_WORLD,MPI_STATUS_IGNORE);

  MPI_Sendrecv(&f2[field_index(0,0,rank_ny-1,0)],transfer_doubles,MPI_DOUBLE,rankp1,rank,&f2[field_index(0,0,-1,0)],
  transfer_doubles,MPI_DOUBLE,rankm1,rankm1,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
  
  MPI_Sendrecv(&f2[field_index(0,0,0,0)],transfer_doubles,MPI_DOUBLE,rankm1,rank,&f2[field_index(0,0,rank_ny,0)],
  transfer_doubles,MPI_DOUBLE,rankp1,rankp1,MPI_COMM_WORLD,MPI_STATUS_IGNORE);

  stream(f1,fc1,f2,fc2,rank_ny);

  compute_rho_u(fc1,fc2,rho1,rho2,ux,uy,uz,rank_ny);

  MPI_Sendrecv(&rho1[sn_index(0,0,rank_ny-1)],transfer_prop,MPI_DOUBLE,rankp1,rank,&rho1[sn_index(0,0,-1)],
  transfer_prop,MPI_DOUBLE,rankm1,rankm1,MPI_COMM_WORLD,MPI_STATUS_IGNORE);

  MPI_Sendrecv(&rho1[sn_index(0,0,0)],transfer_prop,MPI_DOUBLE,rankm1,rank,&rho1[sn_index(0,0,rank_ny)],
  transfer_prop,MPI_DOUBLE,rankp1,rankp1,MPI_COMM_WORLD,MPI_STATUS_IGNORE);

  MPI_Sendrecv(&rho2[sn_index(0,0,rank_ny-1)],transfer_prop,MPI_DOUBLE,rankp1,rank,&rho2[sn_index(0,0,-1)],
  transfer_prop,MPI_DOUBLE,rankm1,rankm1,MPI_COMM_WORLD,MPI_STATUS_IGNORE);

  MPI_Sendrecv(&rho2[sn_index(0,0,0)],transfer_prop,MPI_DOUBLE,rankm1,rank,&rho2[sn_index(0,0,rank_ny)],
  transfer_prop,MPI_DOUBLE,rankp1,rankp1,MPI_COMM_WORLD,MPI_STATUS_IGNORE);

  collide(fc1,fc2,rho1,rho2,ux,uy,uz,solid,rank_ny,rank_ystart);

  double* temp1 = f1;
  f1 = fc1;
  fc1 = temp1;

  double* temp2 = f2;
  f2 = fc2;
  fc2 = temp2;
 }

 rhof = gather(rho1,rank_ny,rank_ystart,size,rank); 

 //visualization//
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
 free(fc1);
 free(fc2);
 free(rho1);
 free(rho2);
 free(ux);
 free(uy);
 free(uz);
 free(rhof);
 free(solid);

 MPI_Finalize();
 return 0;
}