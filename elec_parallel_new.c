#include<stdio.h>
#include<mpi.h>
#include<math.h>
#include<stdlib.h>
#include<string>
#include<sstream>
#include<fstream>

const int NX = 426;
const int NY = 401;
const int ndir = 9;

const double w0 = 4.0/9.0;
const double ws = 1.0/9.0;
const double wd = 1.0/36.0;

const int oxygen_inlet1 = 191;
const int oxygen_inlet2 = 210;

const int water_inlet1 = 400;
const int water_inlet2 = NX;

const int outlet1 = 400;
const int outlet2 = NX;

const double wi[] = {w0,ws,ws,ws,ws,wd,wd,wd,wd};
const int dirx[] = {0,1,0,-1,0,1,-1,-1,1};
const int diry[] = {0,0,1,0,-1,1,1,-1,-1};

const int oopf[] = {0,3,4,1,2,7,8,5,6};

const double rho_w = 1;
const double rho_o = 0.001225;
const double rho_ol = 0.001225;

//const double df = -0.00041;
const double G1 = 2.3;
const double G2 = 2.3;
const double Gads1 = 0.01;
const double Gads2 = -0.01;

const double tau1 = 1.0;
const double tau2 = 1.0;
const int tend = 100;
const double omega1 = 1.0/tau1;
const double omega2 = 1.0/tau2;

inline size_t scalar_index(int x,int y)
{
 return NX*y + x;
}

inline size_t sn_index(int x,int y)
{
 return NX*(y+1) + x;
}

inline size_t field_index(int x,int y,int d)
{
 return ndir*(NX*(y+1) + x) + d;
}

void solid_node(int* solid,int* pm,int size,int rank_ny,int rank_ystart)
{

 for(int y=0,k=rank_ystart; y<rank_ny && k<NY; y++,k++)
 {
  for(int x=0; x<NX; x++)
  {
   solid[sn_index(x,y)] = pm[scalar_index(x,k)];
  }
 }

 for(int y=0,k=rank_ystart; y<rank_ny && k<NY; y++,k++)
 {
  solid[sn_index(0,y)] = 1;
  solid[sn_index(NX-1,y)] = 1;
 }

 for(int y=0,k=rank_ystart; y<rank_ny && k<NY; y++,k++)
 {
  for(int x=0; x<NX; x++)
  {
   if(k==0)
   {
   	solid[sn_index(x,y)]=1;
   }

   if(k==NY-1)
   {
    solid[sn_index(x,y)]=1;
   }
  }
 }

 //oxygen inlet

 for(int y=0,k=rank_ystart; y<rank_ny && k<NY; y++,k++)
 {
   if(k>=oxygen_inlet1 && k<=oxygen_inlet2)
   {
    solid[sn_index(0,y)]=0;
   }
 }


 for(int y=0,k=rank_ystart; y<rank_ny && k<NY; y++,k++)
 {
  for(int x=water_inlet1; x<water_inlet2; x++)
  {
   if(k==0)
   {
    solid[sn_index(x,y)]=0;
   }
  }
 }

 for(int y=0,k=rank_ystart; y<rank_ny && k<NY; y++,k++)
 {
  for(int x=outlet1; x<outlet2; x++)
  {
   if(k==NY-1)
   {
    solid[sn_index(x,y)]=0;
   }
  }
 }

}

void rho_init(double* rho1,double* rho2,int size,int rank,int rank_ny,int rank_ystart)
{
 
 for(int y=0; y<rank_ny; y++)
 {
  for(int x=0; x<NX; x++)
  {
   rho1[sn_index(x,y)]=rho_o;
   rho2[sn_index(x,y)]=rho_w;
  }
 }

 for(int y=0,k=rank_ystart; y<rank_ny && k<NY; y++,k++)
 {
  for(int x=0; x<=10; x++)
  {
   if(k>=oxygen_inlet1 && k<=oxygen_inlet2)
   {
	  rho1[sn_index(x,y)]=rho_w;
    rho2[sn_index(x,y)]=rho_o;
   }
  }
 }

}

void pdf_init(double* rho1,double* rho2,double* f1,double* f2,int size,int rank,int rank_ny,int rank_ystart)
{

  for(int y=0,k=rank_ystart; y<rank_ny && k<NY; y++,k++)
  {
   for(int x=0; x<NX; x++)
   {
 	  for (int d=0; d<ndir; d++)
 	  {
 	   f1[field_index(x,y,d)]=wi[d]*rho1[sn_index(x,y)];
 	   f2[field_index(x,y,d)]=wi[d]*rho2[sn_index(x,y)];
 	  }
   }
  }

}

void compute_rho_u(double* fc1,double* fc2,double* rho1,double* rho2,double* ux,double* uy,int* solid,int rank_ny,int rank_ystart)
{
 
 for(int y=0,k=rank_ystart; y<rank_ny && k<NY; y++,k++)
 {
  for(int x=0; x<NX; x++)
  {
   double r1 = 0.0;
   double r2 = 0.0;
   double u = 0.0;
   double v = 0.0;

   for(int d=0; d<ndir; d++)
   {
    r1 += fc1[field_index(x,y,d)];
    r2 += fc2[field_index(x,y,d)];
   }
    
   for(int d=0; d<ndir; d++) 
   {  
    u += dirx[d]*(fc1[field_index(x,y,d)]+fc2[field_index(x,y,d)]);
    v += diry[d]*(fc1[field_index(x,y,d)]+fc2[field_index(x,y,d)]);
   }

    rho1[sn_index(x,y)] = r1;
    rho2[sn_index(x,y)] = r2;

    ux[scalar_index(x,y)] = u/(r1/tau1 + r2/tau2);
    uy[scalar_index(x,y)] = v/(r1/tau1 + r2/tau2);
   }
  }

}

void macro_boundary_condition(double* rho1,double* rho2,double* ux,double* uy,int rank_ny,int rank_ystart)
{

 for(int y=0,k=rank_ystart; y<rank_ny && k<NY; y++,k++)
 {
  if(k >= oxygen_inlet1 && k <= oxygen_inlet2)
  {
   rho1[sn_index(0,y)]=rho_w;
   rho2[sn_index(0,y)]=rho_o;

   ux[scalar_index(0,y)] = 0.01;
   uy[scalar_index(0,y)] = 0.0;
  }
 }

 for(int y=0,k=rank_ystart; y<rank_ny && k<NY; y++,k++)
 {
  if(k==0)
  {
   for(int x=water_inlet1; x<water_inlet2; x++)
   {
    ux[scalar_index(x,y)] = 0.0;
    uy[scalar_index(x,y)] = 0.00001;
   }
  }

  if(k==NY-1)
  {
   for(int x=outlet1; x<outlet2; x++)
   {
   	rho1[sn_index(x,y)] = rho1[sn_index(x,y-1)];
    rho2[sn_index(x,y)] = rho2[sn_index(x,y-1)];

    ux[scalar_index(x,y)] = 0.0;
    uy[scalar_index(x,y)] = uy[scalar_index(x,y-1)];
   }
  }
 }

}

void forcing(double* rho1,double* rho2,int* solid,double* fx1,double* fx2,double* fsx1,double* fsx2,double* fy1,double* fy2,
  double*  fsy1, double* fsy2,int x,int y,int k)
{
 for (int d=0; d<ndir; d++)
 {
  int xm = (NX+x+dirx[d])%NX;
  int ym = (y+diry[d])%NY;

  *fx1 += wi[d]*dirx[d]*rho2[sn_index(xm,ym)];
  *fy1 += wi[d]*diry[d]*rho2[sn_index(xm,ym)];

  *fx2 += wi[d]*dirx[d]*rho1[sn_index(xm,ym)];
  *fy2 += wi[d]*diry[d]*rho1[sn_index(xm,ym)];

  *fsx1 += wi[d]*dirx[d]*solid[sn_index(xm,ym)];
  *fsy1 += wi[d]*diry[d]*solid[sn_index(xm,ym)];

  *fsx2 += wi[d]*dirx[d]*solid[sn_index(xm,ym)];
  *fsy2 += wi[d]*diry[d]*solid[sn_index(xm,ym)];
 }
  
 *fx1 = -G1 * (*fx1) * rho1[sn_index(x,y)];
 *fy1 = -G1 * (*fy1) * rho1[sn_index(x,y)];

 *fx2 = -G2 * (*fx2) * rho2[sn_index(x,y)];
 *fy2 = -G2 * (*fy2) * rho2[sn_index(x,y)];

 *fsx1 = -Gads1 * (*fsx1) * rho1[sn_index(x,y)];
 *fsy1 = -Gads1 * (*fsy1) * rho1[sn_index(x,y)];

 *fsx2 = -Gads2 * (*fsx2) * rho2[sn_index(x,y)];
 *fsy2 = -Gads2 * (*fsy2) * rho2[sn_index(x,y)];

 if(k>=oxygen_inlet1 && k<=oxygen_inlet2 && x==0)
 {
  x = x+1;
  *fx1 = 0.0;
  *fy1 = 0.0;

  *fx2 = 0.0;
  *fy2 = 0.0;

  for (int d=0; d<ndir; d++)
  {
   int xm = (NX+x+dirx[d])%NX;
   int ym = (y+diry[d])%NY;

   *fx1 += wi[d]*dirx[d]*rho2[sn_index(xm,ym)];
   *fy1 += wi[d]*diry[d]*rho2[sn_index(xm,ym)];

   *fx2 += wi[d]*dirx[d]*rho1[sn_index(xm,ym)];
   *fy2 += wi[d]*diry[d]*rho1[sn_index(xm,ym)];
  }

  *fx1 = -G1 * (*fx1) * rho1[sn_index(x,y)];
  *fy1 = -G1 * (*fy1) * rho1[sn_index(x,y)];

  *fx2 = -G2 * (*fx2) * rho2[sn_index(x,y)];
  *fy2 = -G2 * (*fy2) * rho2[sn_index(x,y)];
 }

 if(x>=water_inlet1 && x<water_inlet2 && k==0)
 {
  y = y+1;
  *fx1 = 0.0;
  *fy1 = 0.0;

  *fx2 = 0.0;
  *fy2 = 0.0;

   for (int d=0; d<ndir; d++)
   {
    int xm = (NX+x+dirx[d])%NX;
    int ym = (y+diry[d])%NY;

    *fx1 += wi[d]*dirx[d]*rho2[sn_index(xm,ym)];
    *fy1 += wi[d]*diry[d]*rho2[sn_index(xm,ym)];

    *fx2 += wi[d]*dirx[d]*rho1[sn_index(xm,ym)];
    *fy2 += wi[d]*diry[d]*rho1[sn_index(xm,ym)];
   }

   *fx1 = -G1 * (*fx1) * rho1[sn_index(x,y)];
   *fy1 = -G1 * (*fy1) * rho1[sn_index(x,y)];

   *fx2 = -G2 * (*fx2) * rho2[sn_index(x,y)];
   *fy2 = -G2 * (*fy2) * rho2[sn_index(x,y)]; 
 }

 if(x>=outlet1 && x<outlet2 && k==NY-1)
 {
   *fx1 = 0.0;
   *fy1 = 0.0;

   *fx2 = 0.0;
   *fy2 = 0.0;
 }
}

void collide(double* f1,double* f2,double* rho1,double* rho2,double* ux,double* uy,int *solid,int rank_ny,int rank_ystart)
{

 double om1 = 1.0 - omega1;
 double om2 = 1.0 - omega2;

 for(int y=0,k=rank_ystart; y<rank_ny && k<NY; y++,k++)
 {
  for(int x=0; x<NX; x++)
  {
   double r1 = rho1[sn_index(x,y)];
   double r2 = rho2[sn_index(x,y)];

   double u = ux[scalar_index(x,y)];
   double v = uy[scalar_index(x,y)];

   double fx1 = 0.0;
   double fx2 = 0.0; 
   double fsx1 = 0.0; 
   double fsx2 = 0.0; 
   double fy1 = 0.0; 
   double fy2 = 0.0; 
   double fsy1 = 0.0;
   double fsy2 = 0.0;

   forcing(rho1,rho2,solid,&fx1,&fx2,&fsx1,&fsx2,&fy1,&fy2,&fsy1,&fsy2,x,y,k);
   
   double ueqx1 = u + (fx1 + fsx1)*tau1/r1;
   double veqy1 = v + (fy1 + fsy1)*tau1/r1;

   double ueqx2 = u + (fx2 + fsx2)*tau2/r2;
   double veqy2 = v + (fy2 + fsy2)*tau2/r2;
   
   if(solid[sn_index(x,y)] == 1)
   {
    double aa[9] , bb[9];

    for(int d=0; d<ndir; d++)
    {
     aa[d] = f1[field_index(x,y,oopf[d])];
     bb[d] = f2[field_index(x,y,oopf[d])];
    }

    for(int d=0; d<ndir; d++)
    {
     f1[field_index(x,y,d)] = aa[d];
     f2[field_index(x,y,d)] = bb[d];
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

     cidotu1 = dirx[d]*ueqx1 + diry[d]*veqy1;

     cidotu2 = dirx[d]*ueqx2 + diry[d]*veqy2;

     feq1 = wi[d]*r1*(1.0 + 3.0*cidotu1 + 4.5*cidotu1*cidotu1 - 1.5*(ueqx1*ueqx1+veqy1*veqy1));

     feq2 = wi[d]*r2*(1.0 + 3.0*cidotu2 + 4.5*cidotu2*cidotu2 - 1.5*(ueqx2*ueqx2+veqy2*veqy2));

     f1[field_index(x,y,d)] = om1*f1[field_index(x,y,d)] + omega1*feq1;
     f2[field_index(x,y,d)] = om2*f2[field_index(x,y,d)] + omega2*feq2;
    } 
   }

  }
 }
}

void stream(double* f1,double* fc1,double* f2,double* fc2,int rank_ny)
{

 for(int y=0; y<rank_ny; y++)
 {
  for(int x=0; x<NX; x++)
  {
   for(int d=0; d<ndir; d++)
   {
    int xmd = (NX+x-dirx[d])%NX;
    int ymd = (y-diry[d])%NY;

    fc1[field_index(xmd,ymd,d)] = f1[field_index(x,y,d)];

    fc2[field_index(xmd,ymd,d)] = f2[field_index(x,y,d)];
   }
  }
 }

}

int* gather_i(int *rho,int num,int rank_ny,int rank_ystart,int size,int rank)
{
 int tag;
 MPI_Status status;
 int *r = NULL;
 int *rank_ny0 = NULL,*rank_ystart0 = NULL;

 if(rank==num)
 {
  r = (int*) malloc(NX*NY*sizeof(int));
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
  MPI_Send(&rho[sn_index(0,0)],NX*rank_ny,MPI_INT,num,tag,MPI_COMM_WORLD);
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
    MPI_Recv(&r[NX*rank_ystart0[i]],NX*rank_ny0[i],MPI_INT,i,tag,MPI_COMM_WORLD,&status);
   }
  }
  
  
  for(int y=0; y<rank_ny; y++)
  {
   for(int x=0; x<NX; x++)
   {
    r[scalar_index(x,y+rank_ystart)] = rho[sn_index(x,y)];
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
  r = (double*) malloc(NX*NY*sizeof(double));
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
  MPI_Send(&rho[sn_index(0,0)],NX*rank_ny,MPI_DOUBLE,num,tag,MPI_COMM_WORLD);
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
    MPI_Recv(&r[NX*rank_ystart0[i]],NX*rank_ny0[i],MPI_DOUBLE,i,tag,MPI_COMM_WORLD,&status);
   }
  }
  
  
  for(int y=0; y<rank_ny; y++)
  {
   for(int x=0; x<NX; x++)
   {
    r[scalar_index(x,y+rank_ystart)] = rho[sn_index(x,y)];
   }
  }
  
 }

 return r;
 
 free(r);
 free(rank_ystart0);
 free(rank_ny0);

}

void SaveTXTImageDataInt(int* rhof,std::string fileName)
{

 FILE* pFile;
 pFile = fopen(fileName.c_str(), "w");

 
 for (int y = 0; y < NY; y++) {
    for (int x = 0; x < NX; x++) {

       fprintf(pFile, "%d\n", rhof[scalar_index(x, y)]);
     }
  }

}

void SaveTXTImageDataDouble(double* rhof,std::string fileName)
{

 FILE* pFile;
 pFile = fopen(fileName.c_str(), "w");

 
 for (int y = 0; y < NY; y++) {
    for (int x = 0; x < NX; x++) {

       fprintf(pFile, "%g\n", rhof[scalar_index(x, y)]);
     }
  }

}

int main(int argc,char *argv[])
{
 int rank,size;
 MPI_Init(&argc,&argv);
 MPI_Comm_size(MPI_COMM_WORLD,&size);
 MPI_Comm_rank(MPI_COMM_WORLD,&rank);

 int rank_ny,rank_ystart;

 if(rank < NY % size)
 {
  rank_ny = NY/size + 1;
  rank_ystart = rank * rank_ny;
 }
 	
 else
 {
  rank_ny = NY/size;
  rank_ystart = NY - (size-rank)*rank_ny;
 }

 int *solid = (int*) calloc(NX*(rank_ny+2),sizeof(int));
 int *solidf = (int*) calloc(NX*NY,sizeof(int));
 int* pm = (int*) calloc(NX*NY,sizeof(int));

 double *f1 = (double*) calloc(NX*(rank_ny+2)*ndir,sizeof(double));
 double *f2 = (double*) calloc(NX*(rank_ny+2)*ndir,sizeof(double));

 double *fc1 = (double*) calloc(NX*(rank_ny+2)*ndir,sizeof(double));
 double *fc2 = (double*) calloc(NX*(rank_ny+2)*ndir,sizeof(double));

 double *rho1 = (double*) calloc(NX*(rank_ny+2),sizeof(double));
 double *rho2 = (double*) calloc(NX*(rank_ny+2),sizeof(double));
 double *rhof = (double*) malloc(NX*NY*sizeof(double));

 double *ux = (double*) calloc(NX*rank_ny,sizeof(double));
 double *uy = (double*) calloc(NX*rank_ny,sizeof(double));

 std::fstream out_p;
 out_p.open("porous.txt");

 for(int y=0; y<NY; ++y)
 {
  for(int x=0; x<NX; ++x)
  { 
   out_p >> pm[scalar_index(x,y)];
  } 
 }
 
 solid_node(solid,pm,size,rank_ny,rank_ystart);

 rho_init(rho1,rho2,size,rank,rank_ny,rank_ystart);

 pdf_init(rho1,rho2,f1,f2,size,rank,rank_ny,rank_ystart);

 int rankp1 = (rank+1) % size;
 int rankm1 = (size+rank-1) % size;

 size_t transfer_doubles = ndir*NX;
 size_t transfer_psi = NX;

 MPI_Sendrecv(&solid[sn_index(0,rank_ny-1)],transfer_psi,MPI_INT,rankp1,rank,&solid[sn_index(0,-1)],
 transfer_psi,MPI_INT,rankm1,rankm1,MPI_COMM_WORLD,MPI_STATUS_IGNORE);

 MPI_Sendrecv(&solid[sn_index(0,0)],transfer_psi,MPI_INT,rankm1,rank,&solid[sn_index(0,rank_ny)],
 transfer_psi,MPI_INT,rankp1,rankp1,MPI_COMM_WORLD,MPI_STATUS_IGNORE);

 for(int t=0; t<tend; t++)
 {
  compute_rho_u(f1,f2,rho1,rho2,ux,uy,solid,rank_ny,rank_ystart);

  MPI_Sendrecv(&rho1[sn_index(0,rank_ny-1)],transfer_psi,MPI_DOUBLE,rankp1,rank,&rho1[sn_index(0,-1)],
  transfer_psi,MPI_DOUBLE,rankm1,rankm1,MPI_COMM_WORLD,MPI_STATUS_IGNORE);

  MPI_Sendrecv(&rho1[sn_index(0,0)],transfer_psi,MPI_DOUBLE,rankm1,rank,&rho1[sn_index(0,rank_ny)],
  transfer_psi,MPI_DOUBLE,rankp1,rankp1,MPI_COMM_WORLD,MPI_STATUS_IGNORE);

  MPI_Sendrecv(&rho2[sn_index(0,rank_ny-1)],transfer_psi,MPI_DOUBLE,rankp1,rank,&rho2[sn_index(0,-1)],
  transfer_psi,MPI_DOUBLE,rankm1,rankm1,MPI_COMM_WORLD,MPI_STATUS_IGNORE);

  MPI_Sendrecv(&rho2[sn_index(0,0)],transfer_psi,MPI_DOUBLE,rankm1,rank,&rho2[sn_index(0,rank_ny)],
  transfer_psi,MPI_DOUBLE,rankp1,rankp1,MPI_COMM_WORLD,MPI_STATUS_IGNORE);

  macro_boundary_condition(rho1,rho2,ux,uy,rank_ny,rank_ystart);

  MPI_Sendrecv(&rho1[sn_index(0,rank_ny-1)],transfer_psi,MPI_DOUBLE,rankp1,rank,&rho1[sn_index(0,-1)],
  transfer_psi,MPI_DOUBLE,rankm1,rankm1,MPI_COMM_WORLD,MPI_STATUS_IGNORE);

  MPI_Sendrecv(&rho1[sn_index(0,0)],transfer_psi,MPI_DOUBLE,rankm1,rank,&rho1[sn_index(0,rank_ny)],
  transfer_psi,MPI_DOUBLE,rankp1,rankp1,MPI_COMM_WORLD,MPI_STATUS_IGNORE);

  MPI_Sendrecv(&rho2[sn_index(0,rank_ny-1)],transfer_psi,MPI_DOUBLE,rankp1,rank,&rho2[sn_index(0,-1)],
  transfer_psi,MPI_DOUBLE,rankm1,rankm1,MPI_COMM_WORLD,MPI_STATUS_IGNORE);

  MPI_Sendrecv(&rho2[sn_index(0,0)],transfer_psi,MPI_DOUBLE,rankm1,rank,&rho2[sn_index(0,rank_ny)],
  transfer_psi,MPI_DOUBLE,rankp1,rankp1,MPI_COMM_WORLD,MPI_STATUS_IGNORE);

  collide(f1,f2,rho1,rho2,ux,uy,solid,rank_ny,rank_ystart);

  MPI_Sendrecv(&f1[field_index(0,rank_ny-1,0)],transfer_doubles,MPI_DOUBLE,rankp1,rank,&f1[field_index(0,-1,0)],
  transfer_doubles,MPI_DOUBLE,rankm1,rankm1,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
  
  MPI_Sendrecv(&f1[field_index(0,0,0)],transfer_doubles,MPI_DOUBLE,rankm1,rank,&f1[field_index(0,rank_ny,0)],
  transfer_doubles,MPI_DOUBLE,rankp1,rankp1,MPI_COMM_WORLD,MPI_STATUS_IGNORE);

  MPI_Sendrecv(&f2[field_index(0,rank_ny-1,0)],transfer_doubles,MPI_DOUBLE,rankp1,rank,&f2[field_index(0,-1,0)],
  transfer_doubles,MPI_DOUBLE,rankm1,rankm1,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
  
  MPI_Sendrecv(&f2[field_index(0,0,0)],transfer_doubles,MPI_DOUBLE,rankm1,rank,&f2[field_index(0,rank_ny,0)],
  transfer_doubles,MPI_DOUBLE,rankp1,rankp1,MPI_COMM_WORLD,MPI_STATUS_IGNORE);

  stream(f1,fc1,f2,fc2,rank_ny);

  double* temp1 = f1;
  f1 = fc1;
  fc1 = temp1;

  double* temp2 = f2;
  f2 = fc2;
  fc2 = temp2;

  if (t%50 == 0 && rank==0)
      printf("%d\n", t);
  
  if(t%10==0)
  { 
   rhof = gather_d(rho2,0,rank_ny,rank_ystart,size,rank);  

   if(rank==0) 
   { 
     std::string densityFileStub("rho");
     std::string fileSuffix1(".txt");
     std::stringstream ts_ind1;
     int vtk_ts1 = t;
     std::string fileName1;
     ts_ind1 << vtk_ts1;

     fileName1 = densityFileStub + ts_ind1.str() + fileSuffix1;
     ts_ind1.str("");
     SaveTXTImageDataDouble(rhof, fileName1); 
   }
  }
 }
 
 MPI_Finalize();
 return 0;

} 