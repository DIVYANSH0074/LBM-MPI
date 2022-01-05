#include<stdio.h>
#include<mpi.h>
#include<math.h>
#include<stdlib.h>
#include<string>
#include<sstream>
#include<fstream>

const int NX = 50;
const int NY = 50;
const int NZ = 60;
const int layer = 40;
const int ndir = 19;

const int radius = 5;
// const int radius2 = 15;
const int start = 15;

const int nx = 0;
const int ny = 0;
const int nz = 0;

const double w0 = 1.0/3.0;
const double ws = 1.0/18.0;
const double wd = 1.0/36.0;

const double wi[] = {w0, ws, ws, ws, ws, ws, ws, wd, wd, wd, wd, wd, wd, wd, wd, wd, wd, wd, wd};

const int dirx[] = {0, 1, -1, 0, 0, 0, 0, 1, 1, 1, 1, -1, -1, -1, -1, 0, 0, 0, 0};
const int diry[] = {0, 0, 0, 1, -1, 0, 0, 1, -1, 0, 0, 1, -1, 0, 0, 1, 1, -1, -1};
const int dirz[] = {0, 0, 0, 0, 0, 1, -1, 0, 0, 1, -1, 0, 0, 1, -1, 1, -1, 1, -1};

const int oopf[] = {0, 2, 1, 4, 3, 6, 5, 12, 11, 14, 13, 8, 7, 10, 9, 18, 17, 16, 15};

const double rho_l = 0.39;
const double rho_v = 0.03;
const double density_air = 7.3160*pow(10,-4);
const double density_solid = 0.50;

const double Tr = (273.15+45.0)/647.0;
const double Tc = 0.09432;
const double T_lu = Tr*Tc;

const double Pc = 0.0044;
const double rho_cr = 0.130;

const double G = -1.0;
const double p_cr = 0.0044;

const double p_g = 1.01325*pow(10,5); // pa
const double p_gstar = 9.595*pow(10,3); // pa
const double kk = (p_g - p_gstar)/(p_g); // ...
const double lx = 4.0*pow(10,-6);
const double delta = 3.316*pow(10,-5); //m2/s
const double mv = 18.0*pow(10,-3);
const double R = 8.314;
const double T = 273.15+45.0;

const double A = lx*lx;
const double Beta = -A*delta*p_g*mv/(R*T*lx);
const double beta = 1.16;
const int dum = 1;
const double blog = Beta*log(kk);

const double tau1 = 1.0;
const double tau2 = 1.0;
const int tend = 500000;
const double omega1 = 1.0/tau1;
const double omega2 = 1.0/tau2;

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

void solid_node(int* solid,int* pm,int size,int rank,int rank_ny,int rank_ystart)
{

 for(int z=0,k=rank_ystart; z<rank_ny && k<10; z++,k++)
 {
  for(int y=0; y<NY; y++)
  {
   for(int x=0; x<NX; x++)
   {
    solid[scalar_index(x,y,z)] = 1;
   }
  }
 }

}

void pdf_init(double* f1,double* f2,int size,int rank,int rank_ny,int rank_ystart)
{
 
 for(int z=0,k=rank_ystart; z<rank_ny && k<NZ; z++,k++)
 {
  for(int y=0; y<NY; y++)
  {
   for(int x=0; x<NX; x++)
   {
    for(int d=0; d<ndir; d++)
    {
     
     if((x-NX/2)*(x-NX/2) + (y-NY/2)*(y-NY/2) + (k-start)*(k-start) <= radius*radius) 
     {
      f1[field_index(x,y,z,d)] = wi[d]*0.39/1.005;
      f2[field_index(x,y,z,d)] = wi[d]*0.0;
     }
     
     else
     {
      f1[field_index(x,y,z,d)] = wi[d]*0.03;
      f2[field_index(x,y,z,d)] = wi[d]*0.0041;
     }

    }
   }
  }
 }

}

void stream(double* f1,double* fc1,double* f2,double* fc2,int rank_ny)
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

void boundary_condition(double* fc1,double* fc2,int rank_ny,int rank_ystart)
{
 
 for(int z=0,k=rank_ystart; z<rank_ny && k<NZ; k++,z++)
 {
  for(int y=0; y<NY; y++)
  {
   for(int x=0; x<NX; x++)
   {
    if(k == NZ-1)
    {
     fc1[field_index(x,y,z,5)] = 0.0;
     fc1[field_index(x,y,z,9)] = 0.0;
     fc1[field_index(x,y,z,13)] = 0.0;
     fc1[field_index(x,y,z,15)] = 0.0;
     fc1[field_index(x,y,z,17)] = 0.0;

     fc2[field_index(x,y,z,5)] = wi[5]*density_air;
     fc2[field_index(x,y,z,9)] = wi[9]*density_air;
     fc2[field_index(x,y,z,13)] = wi[13]*density_air;
     fc2[field_index(x,y,z,15)] = wi[15]*density_air;
     fc2[field_index(x,y,z,17)] = wi[17]*density_air;
    }

    // if(x == 0)
    // {
     
    // }

    // if(x == NX-1)
    // {

    // }

    // if(y == 0)
    // {

    // }

    // if(y == NY-1)
    // {

    // }
   }
  }
 }

}

void compute_rho_u(double* fc1,double* fc2,double* rho1,double* rho2,double* ux,double* uy,double* uz,int* solid,int rank_ny,int rank_ystart)
{
 
 for(int z=0,k=rank_ystart; z<rank_ny && k<NZ; z++,k++)
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
    }

    for(int d=0; d<ndir; d++)
    {

     u += dirx[d]*(fc1[field_index(x,y,z,d)]/tau1 + fc2[field_index(x,y,z,d)]/tau2);
     v += diry[d]*(fc1[field_index(x,y,z,d)]/tau1 + fc2[field_index(x,y,z,d)]/tau2);
     w += dirz[d]*(fc1[field_index(x,y,z,d)]/tau1 + fc2[field_index(x,y,z,d)]/tau2); 

    }

    ux[scalar_index(x,y,z)] = u/(r1/tau1 + r2/tau2);
    uy[scalar_index(x,y,z)] = v/(r1/tau1 + r2/tau2);
    uz[scalar_index(x,y,z)] = w/(r1/tau1 + r2/tau2);

    if(k>0)
    {
     if(solid[scalar_index(x,y,z)] == 1)
     {
      r1 = density_solid;
     }
    }

    rho1[scalar_index(x,y,z)] = r1;
    rho2[scalar_index(x,y,z)] = r2;

   }
  }
 }

}

void compute_psi(double* rho1,double* psi,int rank_ny)
{
 
 for(int z=0; z<rank_ny; z++)
 {
  for(int y=0; y<NY; y++)
  {
   for(int x=0; x<NX; x++)
   {
    double r1 = rho1[scalar_index(x,y,z)];
    double P = r1*T_lu*(1 + r1 + r1*r1 - r1*r1*r1)/pow((1-r1),3) - r1*r1; 
    double P_star = P - r1/3.0;
    psi[sn_index(x,y,z)] = pow(-0.5*P_star,0.5); 
   }
  }
 }

}

void forcing_scheme(double* psi,double* fx1,double* fx2,double* fy1,double* fy2,double* fz1,double* fz2,int x,int y,int z)
{
 
 for(int d=0; d<ndir; d++)
 {
  int xm = (NX+x+dirx[d])%NX;
  int ym = (NY+y+diry[d])%NY;
  int zm = (z+dirz[d])%NZ;

  if(d>0 && d<7)
  {
   *fx1 += 2.0*wi[d]*dirx[d]*psi[sn_index(xm,ym,zm)];
   *fy1 += 2.0*wi[d]*diry[d]*psi[sn_index(xm,ym,zm)];
   *fz1 += 2.0*wi[d]*dirz[d]*psi[sn_index(xm,ym,zm)];  
  }

  else
  {
   *fx1 += wi[d]*dirx[d]*psi[sn_index(xm,ym,zm)];
   *fy1 += wi[d]*diry[d]*psi[sn_index(xm,ym,zm)];
   *fz1 += wi[d]*dirz[d]*psi[sn_index(xm,ym,zm)];
  }
 
  *fx2 += wi[d]*dirx[d]*psi[sn_index(xm,ym,zm)]*psi[sn_index(xm,ym,zm)];
  *fy2 += wi[d]*diry[d]*psi[sn_index(xm,ym,zm)]*psi[sn_index(xm,ym,zm)];
  *fz2 += wi[d]*dirz[d]*psi[sn_index(xm,ym,zm)]*psi[sn_index(xm,ym,zm)];
 }

}

void collide(double *fc1,double *fc2,double* psi,double* rho1,double* rho2,double* ux,double* uy,double* uz,int *solid,int rank_ny,int rank_ystart)
{
 double om1 = 1.0 - omega1;
 double om2 = 1.0 - omega2;

 for(int z=0,k=rank_ystart; z<rank_ny && k<NZ; z++,k++)
 {
  for(int y=0; y<NY; y++)
  {
   for(int x=0; x<NX; x++)
   {
    double r1 = rho1[scalar_index(x,y,z)];
    double r2 = rho2[scalar_index(x,y,z)];

    double u = ux[scalar_index(x,y,z)];
    double v = uy[scalar_index(x,y,z)];
    double w = uz[scalar_index(x,y,z)];

    double fx1 = 0.0;
    double fx2 = 0.0;

    double fy1 = 0.0;
    double fy2 = 0.0;

    double fz1 = 0.0;
    double fz2 = 0.0;

    forcing_scheme(psi,&fx1,&fx2,&fy1,&fy2,&fz1,&fz2,x,y,z);

    fx1 = -G*fx1*psi[sn_index(x,y,z)];
    fy1 = -G*fy1*psi[sn_index(x,y,z)];
    fz1 = -G*fz1*psi[sn_index(x,y,z)];

    fx1 = beta*fx1 - (1.0-beta)*G*fx2/2.0;
    fy1 = beta*fy1 - (1.0-beta)*G*fy2/2.0;
    fz1 = beta*fz1 - (1.0-beta)*G*fz2/2.0;

    double ueq = u + 7.5*fx1/r1;
    double veq = v + 7.5*fy1/r1;
    double weq = w + 7.5*fz1/r1;

    if(solid[scalar_index(x,y,z)] == 1)
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
      
      double cidotu1 = 0.0;
      double cidotu2 = 0.0;

      double feq1 = 0.0;
      double feq2 = 0.0;

      if(k<layer)
      {
       cidotu1 = dirx[d]*ueq + diry[d]*veq + dirz[d]*weq;
       feq1 = wi[d]*r1*(1.0 + 3.0*cidotu1 + 4.5*cidotu1*cidotu1 -1.5*(ueq*ueq+veq*veq+weq*weq));
      }

      else
      {
       feq1 = wi[d]*r1;
      }
      
      cidotu2 = dirx[d]*u + diry[d]*v + dirz[d]*w;
      feq2 = wi[d]*r2*(1.0 + 3.0*cidotu2 + 4.5*cidotu2*cidotu2 -1.5*(u*u+v*v+w*w));

      fc1[field_index(x,y,z,d)] = om1*fc1[field_index(x,y,z,d)] + omega1*feq1;
      fc2[field_index(x,y,z,d)] = om2*fc2[field_index(x,y,z,d)] + omega2*feq2;  
     }
    }
   }
  }
 }

}

int* gather_i(int *rho,int rank_ny,int rank_ystart,int size,int rank)
{
 
 int tag;
 MPI_Status status;
 int *r = NULL;
 int *rank_ny0 = NULL,*rank_ystart0 = NULL;

 if(rank==0)
 {
  r = (int*) malloc(NX*NY*NZ*sizeof(int));
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
  MPI_Send(&rho[scalar_index(0,0,0)],NX*NY*rank_ny,MPI_INT,0,tag,MPI_COMM_WORLD);
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
   MPI_Recv(&r[NX*NY*rank_ystart0[i]],NX*NY*rank_ny0[i],MPI_INT,i,tag,MPI_COMM_WORLD,&status);
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
 }

 return r;
 
 free(r);
 free(rank_ystart0);
 free(rank_ny0);

}

double* gather_d(double *rho,int num,int rank_ny,int rank_ystart,int size,int rank)
{
 int tag;
 MPI_Status status;
 double *r = NULL;
 int *rank_ny0 = NULL,*rank_ystart0 = NULL;

 if(rank==num)
 {
  r = (double*) malloc(NX*NY*NZ*sizeof(double));
  rank_ny0 = (int*) malloc(size*sizeof(int));
  rank_ystart0 = (int*) malloc(size*sizeof(int));
 }

 if(rank != num)
 {
  tag = 1;
  MPI_Send(&rank_ny,1,MPI_INT,num,tag,MPI_COMM_WORLD);
  tag = 2;
  MPI_Send(&rank_ystart,1,MPI_INT,num,tag,MPI_COMM_WORLD);
  tag = 3;
  MPI_Send(&rho[scalar_index(0,0,0)],NX*NY*rank_ny,MPI_DOUBLE,num,tag,MPI_COMM_WORLD);
 }
 
 if(rank==num)
 {
  tag = 1;
  for(int i=0; i<size; i++)
  {
   if(i==num)
   {
    continue;
   }
   else
   { 
    MPI_Recv(&rank_ny0[i],1,MPI_INT,i,tag,MPI_COMM_WORLD,&status);
   }
  }

  tag = 2;
  for(int i=0; i<size; i++)
  {
   if(i==num)
   {
    continue;
   }
   else
   { 
    MPI_Recv(&rank_ystart0[i],1,MPI_INT,i,tag,MPI_COMM_WORLD,&status);
   }
  }
 }

 if(rank==num)
 {
  tag = 3;
  for(int i=0; i<size; i++)
  {
   if(i==num)
   {
    continue;
   }
   else
   { 
    MPI_Recv(&r[NX*NY*rank_ystart0[i]],NX*NY*rank_ny0[i],MPI_DOUBLE,i,tag,MPI_COMM_WORLD,&status);
   }
  }
  
  for(int z=0; z<rank_ny; z++)
  {
   for(int y=0; y<NY; y++)
   {
    for(int x=0; x<NX; x++)
    {
     r[scalar_index(x,y,z+rank_ystart)] = rho[scalar_index(x,y,z)];
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

 printf("rank: %d of %d\n",rank,size);

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

 int *solid = (int*) calloc(NX*NY*rank_ny,sizeof(int));
 int *solidf = (int*) calloc(NX*NY*NZ,sizeof(int));
 int* pm = (int*) calloc(nx*ny*nz,sizeof(int));

 double *psi = (double*) malloc(NX*NY*(rank_ny+2)*sizeof(double));
 
 double *f1 = (double*) malloc(NX*NY*(rank_ny+2)*ndir*sizeof(double));
 double *f2 = (double*) malloc(NX*NY*(rank_ny+2)*ndir*sizeof(double));

 double *fc1 = (double*) malloc(NX*NY*(rank_ny+2)*ndir*sizeof(double));
 double *fc2 = (double*) malloc(NX*NY*(rank_ny+2)*ndir*sizeof(double));

 double *rho1 = (double*) malloc(NX*NY*rank_ny*sizeof(double));
 double *rho2 = (double*) malloc(NX*NY*rank_ny*sizeof(double));
 double *rhof = (double*) malloc(NX*NY*NZ*sizeof(double));

 double *ux = (double*) malloc(NX*NY*rank_ny*sizeof(double));
 double *uy = (double*) malloc(NX*NY*rank_ny*sizeof(double));
 double *uz = (double*) malloc(NX*NY*rank_ny*sizeof(double));

 solid_node(solid,pm,size,rank,rank_ny,rank_ystart);

 pdf_init(f1,f2,size,rank,rank_ny,rank_ystart);

 int rankp1 = (rank+1) % size;
 int rankm1 = (size+rank-1) % size;

 size_t transfer_doubles = ndir*NX*NY;
 size_t transfer_psi = NX*NY;

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

  boundary_condition(fc1,fc2,rank_ny,rank_ystart);

  compute_rho_u(fc1,fc2,rho1,rho2,ux,uy,uz,solid,rank_ny,rank_ystart);

  compute_psi(rho1,psi,rank_ny);

  MPI_Sendrecv(&psi[sn_index(0,0,rank_ny-1)],transfer_psi,MPI_DOUBLE,rankp1,rank,&psi[sn_index(0,0,-1)],
  transfer_psi,MPI_DOUBLE,rankm1,rankm1,MPI_COMM_WORLD,MPI_STATUS_IGNORE);

  MPI_Sendrecv(&psi[sn_index(0,0,0)],transfer_psi,MPI_DOUBLE,rankm1,rank,&psi[sn_index(0,0,rank_ny)],
  transfer_psi,MPI_DOUBLE,rankp1,rankp1,MPI_COMM_WORLD,MPI_STATUS_IGNORE);

  collide(fc1,fc2,psi,rho1,rho2,ux,uy,uz,solid,rank_ny,rank_ystart);

  double* temp1 = f1;
  f1 = fc1;
  fc1 = temp1;

  double* temp2 = f2;
  f2 = fc2;
  fc2 = temp2;

  if(t%100 == 0 && rank == 0)
    printf("time step: %d\n",t);

  if(t%500==0)
  {

   rhof = gather_d(rho1,0,rank_ny,rank_ystart,size,rank);
   solidf = gather_i(solid,rank_ny,rank_ystart,size,rank);

   if(rank == 0)  
   { 
    for(int z=0; z<NZ; z++)
    {
     for(int y=0; y<NY; y++)
     {
      for(int x=0; x<NX; x++)
      {
       if(solidf[scalar_index(x,y,z)] == 1)
       {
        rhof[scalar_index(x,y,z)] = 0;
       }
      }
     }
    }

    std::string densityFileStub("density");
    std::string fileSuffix(".vtk");
    std::stringstream ts_ind;
    std::string ts_ind_str;
    int l_conv_fact = 1;
    int vtk_ts = t;
    std::string fileName;
    std::string dataName("Density");
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
  }

 }

 // rhof = gather_d(rho1,0,rank_ny,rank_ystart,size,rank);

 // if(rank==0) {  

 //    std::string densityFileStub("rho");
 //    std::string fileSuffix(".vtk");
 //    std::stringstream ts_ind;
 //    std::string ts_ind_str;
 //    int l_conv_fact = 1;
 //    int vtk_ts = 0;
 //    std::string fileName;
 //    std::string dataName("Density");
 //    int dims[3];
 //    dims[0] = NX; dims[1] = NY; dims[2] = NZ;
 //    float origin[3];
 //    origin[0] = 0.; origin[1] = 0.; origin[2] = 0.;
 //    int spacing[3];
 //    spacing[0] = l_conv_fact; spacing[1] = l_conv_fact; spacing[2] = l_conv_fact;

 //   //Writing VTK File
 //    ts_ind << vtk_ts;

 //    fileName = densityFileStub + ts_ind.str() + fileSuffix;
 //    ts_ind.str("");
 //    SaveVTKImageData_ascii(rhof, fileName, dataName, origin, spacing, dims); }

 free(f1);
 free(f2);
 free(fc1);
 free(fc2);
 free(rho1);
 free(rho2);
 free(rhof);
 free(ux);
 free(uy);
 free(uz);
 free(solid);
 free(solidf);
 free(psi);

 MPI_Finalize();
 return 0;

}