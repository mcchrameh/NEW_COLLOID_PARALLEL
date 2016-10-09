#include "CDS_Colloids.h"
#include <random>
#include <time.h>


#define Random_min  -0.05
#define Random_max  0.05 








Colloid::Colloid():BASE()
        {
          Posx.resize(Nparticles);
	  Posy.resize(Nparticles);
	  Fx.resize(Nparticles);
	  Fy.resize(Nparticles);
	  Vx.resize(Nparticles);
	  Vy.resize(Nparticles);



	}
Colloid::~Colloid()
{

 printf("Destructor of Colloid class called\n");

}
void Colloid::Exchange_BcColloids(MPI_Comm new_comm)
{
    MPI_Comm_rank(new_comm, &my2drank);
    MPI_Comm_size(new_comm, &size);



 





}
void Colloid::ComputeForce(MPI_Comm new_comm)//, BODY* bd, Node* root, double dist=0.0)
{
 MPI_Comm_rank(new_comm, &my2drank);
 MPI_Comm_size(new_comm, &size);

 double G4(0.0), AA(0.0), BB(0.0), NEW1(0.0), NEW2(0.0), G2(0.0), G3(0.0), G3A(0.0),C(0.0);
 double DSUMX(0.0),DSUMY(0.0), DSUMXY(0.0);
 BASE::SetEffectiveRadius( new_comm);

   // Coupling( new_comm);

     
   //printf("rank=%d, END INDEX  IN EACH CORE=%d\n",my2drank,nlocal_particles+count_particles);
 for(int i=N_start;i<N_start+ nlocal_particles-1;i++)
    {
     for(int j=i+1;j<N_start+ nlocal_particles;j++)
        {
    //	printf("Colloid (dx=%lf, dy=%lf) \n",bd[i].x, bd[j].x);       
         double dx=bd[i].r[0]-bd[j].r[0];
         double dy=bd[i].r[1]-bd[j].r[1];
      	//printf("Colloid (dx=%lf, dy=%lf) \n",dx, dy);
         //dx=dx-nlocalx*(int((dx/(nlocalx/2)+3)/2)-1.0);
         dy=dy-nlocaly*(int((dy/(nlocaly/2)+3)/2)-1.0);
         double R2=sqrt(pow(dx,2)+pow(dy,2));
         if(R2<R12[i][j])
           {
            Dx[i][j]=dx;
            Dy[i][j]=dy;
	    RI[i][j]=R2;
	    RI[j][i]=R2;
	    C = U1/R12[i][j];
	    AA=RI[i][j]/R12[i][j];
	    G2 = pow(AA + beta,2.0);
	    G3 = 1.0 + ALPHA*(AA + beta);
	    BB = (AA-1.0)*(-ALPHA);
	    G3A = exp(BB);
	    NEW1=pow((1.0 + (ENNE/ALPHA)+ beta),2.0);
	    NEW2=(exp(-ENNE))*((1.0 + ALPHA*(beta + (ENNE/ALPHA)+1.0)))/(NEW1);
	    G4 = ((C*(G3A)*G3)/G2)-(C*NEW2);
            Fx[i]=Fx[i] + G4*(Dx[i][j]/RI[i][j]);
            Fy[i]=Fy[i] + G4*(Dy[i][j]/RI[i][j]);
	//   printf("Colloid Force (fx=%lf, fy=%lf) \n", Fx[i], Fy[i]);
	   // Fx[j]=Fx[j] + G4*(Dx[j][i])/RI[i][j];
	  // Fy[j]=Fy[j] + G4*(Dy[j][i])/RI[i][j];
	    }
        }
     }
//Compute random force
    double csi1(0.0), csi2(0.0);
    srand(time(NULL) + rank);
    double GAMM(0.0);

   for(int i=N_start;i<N_start + nlocal_particles; i++)
      {
       double range =Random_max-Random_min;
       double div =RAND_MAX / range;
       csi1 = Random_min + (rand()/div);
       csi2 = Random_min + (rand()/div);
     //  printf("Random number (csi1=%lf, csi2=%lf) \n", csi1, csi2);
      if(bd[i].index==0)
        {
  	  GAMM= GAMMs;
          bd[i].v[0]= (1.0/GAMM)*(Fx[i]+ TOTX[i])*delta_t + sqrt(2.0*temp1/GAMM)*csi1*sqrt(delta_t);
          bd[i].v[1]= (1.0/GAMM)*(Fy[i]+ TOTY[i])*delta_t + sqrt(2.0*temp1/GAMM)*csi1*sqrt(delta_t);
     	 // printf("Colloid Force (fx=%lf, fy=%lf) \n", Fx[i], Fy[i]);
        }
     if(bd[i].index==1)
	{
         GAMM= GAMMb;
         bd[i].v[0]= (1.0/GAMM)*(Fx[i]+ TOTX[i])*delta_t + sqrt(2.0*temp1/GAMM)*csi2*sqrt(delta_t);
         bd[i].v[1]= (1.0/GAMM)*(Fy[i]+ TOTY[i])*delta_t + sqrt(2.0*temp1/GAMM)*csi2*sqrt(delta_t);
        // printf("Colloid (TOTX =%lf, TOTY=%lf) \n", TOTX[i], TOTY[i]);
   	}
       //position update
       bd[i].r[0] =bd[i].r[0] + bd[i].v[0]; //<---bd[i].vx is already multiplied by t. 
       bd[i].r[1] =bd[i].r[1] + bd[i].v[1];
   //  printf("bd[%d].x=%lf, bd[%d].y=%lf\n",i,bd[i].x,i,bd[i].y);
      //bc condition
     //bd[i].x=bd[i].x-nlocalx*((int)(bd[i].x/nlocalx + 1.0)-1.0);
      if((bd[i].r[0]>(double)nlocalx)||(bd[i].r[1]>(double)nlocaly))
	{
         //std::vector<BODY>::iterator it;
         printf("bd[%d].x=%lf, bd[%d].y=%lf\n",i,bd[i].r[0],i,bd[i].r[1]);		
        // Send_bodies.push_back(bd[i]);
  //	  bd.erase(bd.begin()+i-1);
	}
      //bd[i].x=bd[i].x-(nlocalx/2)*((double)rand())/(RAND_MAX);

      //printf("Rank=%d,Positionx=%lf\n ", my2drank, bd[i].x);
     //bd[i].y=bd[i].y-nlocaly*((int)(bd[i].y/nlocaly + 1.0)-1.0);
     //computing mean displacement
     Posx[i]=Posx[i] + bd[i].v[0];
     Posy[i]=Posy[i] + bd[i].v[1];
     DSUMX  = DSUMX + pow(Posx[i],2.0);
     DSUMY  = DSUMY + pow(Posy[i],2.0);
     DSUMXY = DSUMXY + Posx[i]*Posy[i];
    // printf("Colloid positions (x=%lf, y=%lf) \n", bd[i].x, bd[i].y);
   }
   // DSUMX=DSUMX/
}
/*

void Colloid::WriteToFile_colloid(MPI_Comm new_comm)
{
   MPI_Status status;
   MPI_File     fh;
   MPI_Datatype filetype;
   int dims[2],coords[2], periods[2], start_indices[2], localsizes[2];
   int globalsizes[2];
   globalsizes[0]=BASE::Nx;
   globalsizes[1]=BASE::Ny;
   localsizes[0]=(BASE::Nx/BASE::Procx); localsizes[1]=(BASE::Ny/BASE::Procy) ;
   int local_particle_size =BASE::nlocal_particles;
  
   MPI_Cart_get (new_comm,2,dims,periods, coords );
   MPI_Comm_rank(new_comm, &my2drank);
   MPI_Comm_size(new_comm, &size);
   MPI_Cart_coords(new_comm,my2drank,2,coords);
   start_indices[0]=coords[0]*localsizes[0];
   start_indices[1]=coords[1]*localsizes[1];
   printf("rank=%d,r=%d,start_indices=(%d,%d)\n",my2drank,Nx%Procx,start_indices[0],start_indices[1]);
   MPI_Type_create_subarray(2, globalsizes, localsizes, start_indices,MPI_ORDER_C, MPI_DOUBLE, &filetype);
   MPI_Type_commit(&filetype);
   
   char outputfilename[]="datafile_colloid";
   char filememview[]="native";
   MPI_File_open(new_comm, outputfilename, MPI_MODE_CREATE | MPI_MODE_WRONLY,MPI_INFO_NULL, &fh);
   MPI_File_set_view(fh, 0, MPI_DOUBLE, filetype, filememview,MPI_INFO_NULL);
   MPI_File_write_all(fh, PHI_local_result,local_array_size , MPI_DOUBLE, &status);
   MPI_File_close(&fh);
   MPI_Type_free(&filetype);


   

}
*/
