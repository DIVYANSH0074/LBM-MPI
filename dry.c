#include<stdio.h>
#include<mpi.h>
#include<math.h>
#include<stdlib.h>
#include<string>
#include<sstream>
#include<fstream>

const int NX = 40;
const int NY = 40;
const int NZ = 300;
const int layer = 200;
const int ndir = 19;

const int nx = 40;
const int ny = 40;
const int nz = 200;

const double w0 = 1.0/3.0;
const double ws = 1.0/18.0;
const double wd = 1.0/36.0;

const double wi[] = {w0, ws, ws, ws, ws, ws, ws, wd, wd, wd, wd, wd, wd, wd, wd, wd, wd, wd, wd};

const int dirx[] = {0, 1, -1, 0, 0, 0, 0, 1, 1, 1, 1, -1, -1, -1, -1, 0, 0, 0, 0};
const int diry[] = {0, 0, 0, 1, -1, 0, 0, 1, -1, 0, 0, 1, -1, 0, 0, 1, 1, -1, -1};
const int dirz[] = {0, 0, 0, 0, 0, 1, -1, 0, 0, 1, -1, 0, 0, 1, -1, 1, -1, 1, -1};

const int oopf[] = {0,2,1,4,3,6,5,12,11,14,13,8,7,10,9,18,17,16,15};

const double rho_l = 0.39;
const double rho_v = 0.03;
const double density_air = 7.3160*pow(10,-4);

const double Tr = 0.4915;
const double Tc = 0.0943;
const double T_lu = Tr*Tc;

const double Pc = 0.0044;
const double rho_cr = 0.130;

//const double G = -1;
const double p_cr = 0.0044;

const double p_g = 1.013*pow(10,5);
const double lx = NX*4.0*pow(10,-6);
const double delta = 0.3205;
const double mv = 18.0*pow(10,-3);
const double R = 8.314;
const double T = 333.15;

const double A = NY*4.0*pow(10,-6)*lx;
const double Beta = -A*delta*p_g*mv/(R*T);
const double beta = 1.16;
const int dum = 1;

const double tau1 = 1.0;
const double tau2 = 1.0;
const int tend = 3;
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

inline size_t sign(double P)
{
 if(P>0.0)
 {
  return 1;
 }
 
 else if(P<0.0)
 {
  return -1;
 }

 else
 {
  return 0;
 }
}

inline double absolute(double P)
{
 if(P<0.00)
 {
  return P*-1.0;
 }

 else
 {
  return P;
 }
}

void solid_node(int* solid,int* pm,int size,int rank,int rank_ny,int rank_ystart)
{
 
 for(int z=0,k=rank_ystart; z<rank_ny && k<layer; z++,k++)
 {
  for(int y=0; y<NY; y++)
  {
   for(int x=0; x<NX; x++)
   {
    solid[scalar_index(x,y,z)] = pm[scalar_index(x,y,k)];
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
     
     if(k<layer) 
     {
      f2[field_index(x,y,z,d)] = wi[d]*0.0;
     }
     
     else
     {
      f2[field_index(x,y,z,d)] = wi[d]*0.0041;
     }

    }
   }
  }
 }

 for(int z=0,k=rank_ystart; z<rank_ny && k<NZ; z++,k++)
 {
  for(int y=0; y<NY; y++)
  {
   for(int x=0; x<NX; x++)
   {
    for(int d=0; d<ndir; d++)
    {
     
     if(k>0 && k<layer) 
     {
      f1[field_index(x,y,z,d)] = wi[d]*0.39/1.005;
     }

     else
     {
      f1[field_index(x,y,z,d)] = wi[d]*0.03;
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
     if(solid[scalar_index(x,y,z)] == 1 || solid[scalar_index(x,y,z)] == 2)
     {
      r1 = 0.39;
     }
    }

    rho1[scalar_index(x,y,z)] = r1;
    rho2[scalar_index(x,y,z)] = r2;

   }
  }
 }

}

void compute_psi(double* rho1,double* psi,int* G_m,int rank_ny)
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
    G_m[scalar_index(x,y,z)] = sign(P_star);

    double P_abs = absolute(P_star); 
    psi[sn_index(x,y,z)] = pow(2.0*P_abs,0.5); 
   }
  }
 }

}

void forcing_scheme(double* psi,double* fx1,double* fx2,double* fy1,double* fy2,double* fz1,double* fz2,int x,int y,int z,int rank_ny)
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

void collide(double *fc1,double *fc2,double* psi,int* G_m,double* rho1,double* rho2,double* ux,double* uy,double* uz,int *solid,int rank_ny,int rank_ystart)
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
    int G = G_m[scalar_index(x,y,z)];

    double u = ux[scalar_index(x,y,z)];
    double v = uy[scalar_index(x,y,z)];
    double w = uz[scalar_index(x,y,z)];

    double fx1 = 0.0;
    double fx2 = 0.0;

    double fy1 = 0.0;
    double fy2 = 0.0;

    double fz1 = 0.0;
    double fz2 = 0.0;

    forcing_scheme(psi,&fx1,&fx2,&fy1,&fy2,&fz1,&fz2,x,y,z,rank_ny);

    fx1 = -G*fx1*psi[sn_index(x,y,z)];
    fy1 = -G*fy1*psi[sn_index(x,y,z)];
    fz1 = -G*fz1*psi[sn_index(x,y,z)];

    fx1 = beta*fx1 - (1.0-beta)*G*fx2/2.0;
    fy1 = beta*fy1 - (1.0-beta)*G*fy2/2.0;
    fz1 = beta*fz1 - (1.0-beta)*G*fz2/2.0;

    double ueq = u + 2.0*fx1/r1;
    double veq = v + 2.0*fy1/r1;
    double weq = w + 2.0*fz1/r1;

    if(solid[scalar_index(x,y,z)] == 1 || solid[scalar_index(x,y,z)] == 2)
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

        fprintf(pFile, " %f  ", rho1[scalar_index(x, y, z)]);
      }
    }
    fprintf(pFile, "\n");
  }

  fclose(pFile);
}

void SaveTXTImageData(double* rhof,std::string fileName)
{

 FILE* pFile;
 pFile = fopen(fileName.c_str(), "w");

 for (int z = 0; z < NZ; z++) {
    for (int y = 0; y < NY; y++) {
      for (int x = 0; x < NX; x++) {

        fprintf(pFile, "%g\n", rhof[scalar_index(x, y, z)]);
      }
    }
  }

}

void SaveTXTImageData_i(int* rhof,std::string fileName)
{

 FILE* pFile;
 pFile = fopen(fileName.c_str(), "w");

 for (int z = 0; z < NZ; z++) {
    for (int y = 0; y < NY; y++) {
      for (int x = 0; x < NX; x++) {

        fprintf(pFile, "%d\n", rhof[scalar_index(x, y, z)]);
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

void total(double* total_water,double* length,double* rho1,int* solid,int rank_ny,int rank_ystart,int rank)
{
 
 for(int z=0,k=rank_ystart; z<rank_ny && k<NZ; z++,k++)
 {
  for(int y=0; y<NY; y++)
  {
   for(int x=0; x<NX; x++)
   {
    double rho_3 = rho1[scalar_index(x,y,z)]*(1.0 - solid[scalar_index(x,y,z)]);

    if(k<layer)
    {
     if(rho_3 > 0.3)
     {
      *total_water = *total_water + 1.0;
     }
    }
    
    if(solid[scalar_index(x,y,z)] == 1 || solid[scalar_index(x,y,z)] == 2)
    {
     *length += 1.0;
    }
    
   }
  }
 }

}

void saturation(double* gtotal_water,double* glength,double* sat,int size)
{
 double tl =0.0;
 double tw = 0.0;

 for(int i=0; i<size; i++)
 {
  tw += gtotal_water[i];
  tl += glength[i]; 
 }

 double voidspace = (double)NX*NY*layer - tl;
 *sat = tw*100.0/voidspace;
  printf("%g\n",*sat);
}

void MV_init(double* rdr,double* rho1,int* G_m,double* psi,int rank_ny,int rank_ystart)
{
 double p_water_c = 31.9;

 for(int z=0,k=rank_ystart; z<rank_ny && k<NZ; z++,k++)
  {
   for(int y=0; y<NY; y++)
   {
    for(int x=0; x<NX; x++)
    {
     rdr[scalar_index(x,y,z)] = (1.0/3.0)*rho1[scalar_index(x,y,z)]*220*pow(10,5)/(p_water_c);
    }
   }
  }

}

void MV(double* rdr,double* Mv,int rank_ny,int rank_ystart)
{
 double l = 0.0;

 for(int z=0,k=rank_ystart; z<rank_ny && k<NZ; z++,k++)
 {
  for(int y=0; y<NY; y++)
  {
   for(int x=0; x<NX; x++)
   {
    if(k == NZ-3)
    { 
     double kk = (p_g - rdr[scalar_index(x,y,z)])/(p_g - rdr[scalar_index(x,y,z+1)]);

     l += Beta*log(kk);
    }
   }
  }
 }

 *Mv = l*3600;

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
 double sat = 0.0, length = 0.0, total_water = 0.0, Mv = 0.0;
 double *glength, *gtotal_water;
 int* G_m = (int*) malloc(NX*NY*rank_ny*sizeof(int));
 double *psi = (double*) malloc(NX*NY*(rank_ny+2)*sizeof(double));
 
 double *f1 = (double*) malloc(NX*NY*(rank_ny+2)*ndir*sizeof(double));
 double *f2 = (double*) malloc(NX*NY*(rank_ny+2)*ndir*sizeof(double));

 double *fc1 = (double*) malloc(NX*NY*(rank_ny+2)*ndir*sizeof(double));
 double *fc2 = (double*) malloc(NX*NY*(rank_ny+2)*ndir*sizeof(double));

 double *rho1 = (double*) malloc(NX*NY*rank_ny*sizeof(double));
 double *rho2 = (double*) malloc(NX*NY*rank_ny*sizeof(double));
 double *rdr = (double*) malloc(NX*NY*rank_ny*sizeof(double));
 double *rhof = (double*) malloc(NX*NY*NZ*sizeof(double));

 double *ux = (double*) malloc(NX*NY*rank_ny*sizeof(double));
 double *uy = (double*) malloc(NX*NY*rank_ny*sizeof(double));
 double *uz = (double*) malloc(NX*NY*rank_ny*sizeof(double));

 if(rank==1)
 {
  glength = (double*) malloc(size*sizeof(double));
  gtotal_water = (double*) malloc(size*sizeof(double));
 }

 std::fstream out_p;
 out_p.open("porous.txt");

 for(int z=0; z<nz; z++)
 {
  for(int y=0; y<ny; ++y)
  {
   for(int x=0; x<nx; ++x)
   { 
    out_p >> pm[scalar_index(x,y,z)];
   } 
  } 
 }

 solid_node(solid,pm,size,rank,rank_ny,rank_ystart);

 pdf_init(f1,f2,size,rank,rank_ny,rank_ystart);

 int rankp1 = (rank+1) % size;
 int rankm1 = (size+rank-1) % size;

 size_t transfer_doubles = ndir*NX*NY;
 size_t transfer_psi = NX*NY;

 std::ofstream out_Mv,out_sat;
 out_Mv.open("Mv_o.txt", std::ios::trunc);
 out_sat.open("sat.txt", std::ios::trunc);

 for(int t=0; t<tend; t++)
 {

  sat = 0.0, length = 0.0, total_water = 0.0, Mv = 0.0;

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

  compute_psi(rho1,psi,G_m,rank_ny);

  MPI_Sendrecv(&psi[sn_index(0,0,rank_ny-1)],transfer_psi,MPI_DOUBLE,rankp1,rank,&psi[sn_index(0,0,-1)],
  transfer_psi,MPI_DOUBLE,rankm1,rankm1,MPI_COMM_WORLD,MPI_STATUS_IGNORE);

  MPI_Sendrecv(&psi[sn_index(0,0,0)],transfer_psi,MPI_DOUBLE,rankm1,rank,&psi[sn_index(0,0,rank_ny)],
  transfer_psi,MPI_DOUBLE,rankp1,rankp1,MPI_COMM_WORLD,MPI_STATUS_IGNORE);

  collide(fc1,fc2,psi,G_m,rho1,rho2,ux,uy,uz,solid,rank_ny,rank_ystart);

  total(&total_water,&length,rho1,solid,rank_ny,rank_ystart,rank);

  MPI_Gather(&total_water,1,MPI_DOUBLE,gtotal_water,1,MPI_DOUBLE,1,MPI_COMM_WORLD);

  MPI_Gather(&length,1,MPI_DOUBLE,glength,1,MPI_DOUBLE,1,MPI_COMM_WORLD);

  if(rank==1)
  {
   /***************** saturation function *********************/

   saturation(gtotal_water,glength,&sat,size);

   std::ofstream out_sat;
   
   out_sat.open("sat.txt", std::ios::app);
 
   out_sat<< sat << std::endl;

   out_sat.close();
  }

  MV_init(rdr,rho1,G_m,psi,rank_ny,rank_ystart);

  MV(rdr,&Mv,rank_ny,rank_ystart);
   
  if(rank == size-1)
  {
   std::ofstream out_Mv;
   out_Mv.open("Mv_o.txt", std::ios::app);
   out_Mv << Mv << std::endl;
   out_Mv.close();
  }

  MPI_Bcast(&sat,1,MPI_DOUBLE,1,MPI_COMM_WORLD);

  if(sat <= 0.0)
  {
   break;
  }

  double* temp1 = f1;
  f1 = fc1;
  fc1 = temp1;

  double* temp2 = f2;
  f2 = fc2;
  fc2 = temp2;

  if(t%100 == 0 && rank == 0)
    printf("time step: %d\n",t);

  if(t%1000==0)
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

       if(solidf[scalar_index(x,y,z)] == 2)
       {
        rhof[scalar_index(x,y,z)] = 0;
       }

       if(solidf[scalar_index(x,y,z)] == 1)
       {
        rhof[scalar_index(x,y,z)] = 0;
       }
      }
     }
    }

    std::string densityFileStub("density_pln");
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

  if(t%50000 == 0)
  {
   rhof = gather_d(rho1,3,rank_ny,rank_ystart,size,rank);

   if(rank == 3)
   {
    std::string densityFileStub("rho1");
    std::string fileSuffix1(".txt");
    std::stringstream ts_ind1;
    int vtk_ts1 = t;
    std::string fileName1;
    ts_ind1 << vtk_ts1;

    fileName1 = densityFileStub + ts_ind1.str() + fileSuffix1;
    ts_ind1.str("");
    SaveTXTImageData(rhof, fileName1);
   }

  
   }

  if(t%50000 == 0)
  {
   rhof = gather_d(rho2,2,rank_ny,rank_ystart,size,rank);

   if(rank == 2)
   {
     std::string densityFileStub("rho2");
     std::string fileSuffix1(".txt");
     std::stringstream ts_ind1;
     int vtk_ts1 = t;
     std::string fileName1;
     ts_ind1 << vtk_ts1;

     fileName1 = densityFileStub + ts_ind1.str() + fileSuffix1;
     ts_ind1.str("");
     SaveTXTImageData(rhof, fileName1);
   }
  }
 } 

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
 free(G_m);
 free(psi);

 MPI_Finalize();
 return 0;
}