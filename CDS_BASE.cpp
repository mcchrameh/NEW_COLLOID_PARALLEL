#include "CDS_BASE.h"
#include <iostream>
#include <random>
#include <stdio.h>
#include <stdlib.h>
#include <cstring>
#include <fstream>
#include <sstream>
#include <string>
#include <time.h>

//#include "mpi.h"
#define MAX_Partices  100000
#define Random_min  -0.05
#define Random_max  0.05 
#define MAX_LINE_LENGTH 200
#define gamma(i,j)     gamma[i][j]
#define PHI_old(i,j)   PHI_old[i][j]
#define PHI(i,j)       PHI[i][j]
#define Laplacian2(i,j) Laplacian2[i][j]
#define POLYMER(i,j)    POLYMER[i*Ny+j]
#define PP3(i,j)       PP3[i*Ny + j]
#define PPP(i,j)       PPP[i*Ny + j]

#define P2(i,j)    P2[i*Ny +j]
#define CP4(i,j)  CP4[i][j]
BASE::BASE()
     {
      MPI_Comm_rank(MPI_COMM_WORLD,&rank);
      MPI_Comm_size(MPI_COMM_WORLD, &size);      
       double rationals[26];
       int integers[8];
       if(rank==0)
	  {
	    Read_input_parameters(integers, rationals);
	   }
       MPI_Bcast(integers,8,MPI_INT,0,MPI_COMM_WORLD);
       MPI_Bcast(rationals,26,MPI_DOUBLE,0,MPI_COMM_WORLD);
       Nx   = integers[0];     Ny    = integers[1]; //Nz=integers[2];
       Procx=integers[2];      Procy = integers[3]; //Procz =integers[5];
       Nbig =integers[4];      Nsmall= integers[5];
       Nparticles=integers[6];
       M = rationals[0]; u=rationals[1]; B = rationals[2];
       D = rationals[3]; A=rationals[4]; v = rationals[5];
       tau= rationals[6]; F=rationals[7]; r=rationals[8];
       Max_time_iteration=rationals[9];
       U1=rationals[10]; CSI=rationals[11];
       ENNE=rationals[12];      ALPHA =rationals[13];
       beta =rationals[14];     GAMMs =rationals[15];
       GAMMb =rationals[16];    temp1 =rationals[17];
       R1small=rationals[18];   R1big =rationals[19];
       ALP1 =rationals[20];     dxx   =rationals[21];
       EMME =rationals[22];     sigma =rationals[23];
       P00s =rationals[24];     P00b  =rationals[25];

        delta_t=0.01; 
       //Nbig = integers[4];
       //Nparticles=integers[5];
       RR1small= R1small*pow((1.0 +(1.0/log(2.0))),(1.0/ALP1));
       RR1big=R1big*pow((1.0 +(1.0/log(2.0))),(1.0/ALP1));
       RCUT = RR1big*0.01; //<---This is arbitrary
       //printf("RCUT=%lf\n",RCUT);

      int coords[2], dims[2], periods[2];

      MPI_Comm new_comm =CreatCartesianTopology();
      MPI_Comm_rank(new_comm,&my2drank);
      MPI_Cart_get(new_comm,2,dims,periods,coords);
      assert(Procx*Procy==size);
      printf("dimensions of topology=(%d,%d)\n",dims[0],dims[1]);
      
      if(Nx%Procx!=0)
	{
         if(coords[0]<(dims[0]-1))
           {
             nlocalx=(int)Nx/Procx;
	   }
	 else
          {
            nlocalx = (int)Nx/Procx + (int)Nx%Procx;

          }
	}
      else
	{
		
          nlocalx = (int)Nx/Procx;
        }

      if(Ny%Procy!=0)
	{
	  if(coords[1]<(dims[1]-1))
	    {
	      nlocaly = (int)Ny/Procy;
	    }
	 else
	    {
	      nlocaly = (int)Ny/Procy + (int)Ny%Procy;
	    }
	 }
      else
	 {
             nlocaly = (int)Ny/Procy;

          }
   
//----------------------------Calculate the local grid on a process------------------------------------//
      locala_x = 0.0 + coords[0]*nlocalx*(1.0/Nx); //hx=1.0/Nx, [0.0, 1.0]x[0.0,1.0], global a =0.0, global b =1.0
      //localb_x= 


/*--Compute Number of linked cells in x and y direction -------*/
    /*  LinkedCells.resize(2);
      CellSize.resize(2);
      for (int direction=0; direction<2; direction++) 
          {
	    LinkedCells[direction] = (int) (al[direction]/RCUT) ; 
	    CellSize[direction]    = (double)(al[direction]/LinkedCells[direction]);
	    printf("Cellsize[%d]=%lf\n",direction, LinkedCells[direction]);
          }
      //---Add ghost cells to number of linked cells---/
      for (int direction=0; direction<2; direction++)
	   {
              LinkedCells[direction]=LinkedCells[direction] + 2;
	   }

      head.resize(LinkedCells[0]*LinkedCells[1]);
      LinkedCell_list.resize(Nparticles);
     */
    //-----------------------------------------------------------------
    //Establish the local number of particles on each MPI process
    //---------------------------------------------------------------
           /*
	   nbig_local= (int)(Nbig/size);
	   nbig_local=(my2drank<Nbig%size)?nbig_local +1 : nbig_local ; //<-----check this
	   nsmall_local= (int)(Nsmall/size);
	   nsmall_local=(my2drank<Nsmall%size)?nsmall_local +1 : nsmall_local ;

           nlocal_particles =Nparticles/size;
           N_start = (int)my2drank*nlocal_particles;
	   if(Nparticles!=(Nbig + Nsmall))
	   {
            printf("particles don't match; Nbig=%d, Nsmall=%d, Nparticles=%d\n", Nbig, Nsmall, Nparticles);
	    printf("Condition Nparticles = Nbig + Nsmall, is not satisfied");
	    assert(Nparticles==(Nbig + Nsmall));
	   }
	   */
       blocklengths[0]=1;
       offsets[0]=0;
       types[0]=MPI_INT;
       MPI_Type_extent(MPI_INT, &extent);
       offsets[1]=1*extent;
       types[1]=MPI_DOUBLE;
       blocklengths[1]=9;
       MPI_Type_struct(2, blocklengths, offsets, types, &particletype);
       MPI_Type_commit(&particletype);
           nlocal_particles =Nparticles/size;
           
	   N_start = (int)my2drank*nlocal_particles;
      /* if(Nparticles%size!=0)
	{
	  if(coords[0]<(dims[0]-1))
	    {
	      nlocal_particles = (int)Nparticles/size;
	      N_start = (int)my2drank*nlocal_particles;

	    }
	 else
	    {
	      nlocal_particles = (int)Nparticles/size + (int)Nparticles%size;
	      N_start = (int)my2drank*nlocal_particles;
	      
	    }
	 }
      else
	 {
             nlocal_particles = (int)Nparticles/size;
             N_start = (int)my2drank*nlocal_particles;

          }
         */
  

	   if(Nparticles%size)
	    {
             nlocal_particles=(my2drank<Nparticles%size)?nlocal_particles +1 : nlocal_particles ;
	     N_start =N_start - ((my2drank < ((Nparticles)%size) ) ? my2drank : ((Nparticles)%size));
	    }
	   printf("Rank=%d, N_start=%d, nlocalParticles=%d\n",my2drank, N_start,nlocal_particles);
      //     printf("Rank=%d,Nlocalx=%d, Nlocaly=%d, total=%d\n", my2drank, nlocalx,nlocaly,nlocalx*nlocaly);	    
            

  
  recvcounts=new int[size];
  displacement=new int[size];
  recvcounts_polymer=new int[size];
  disp_polymer      =new int[size];
 
  disp_polymer[0]=0;
  //recvcounts_polymer[0]=nlocalx*nlocaly;//=(Nx/Procx)*(Ny/Procy);//nlocalx*nlocaly;
  for(int i=0;i<size;i++)
     {
     recvcounts_polymer[i]=((i<(size-1))?(int)Nx/Procx:(int)Nx/Procx + (int)Nx%Procx)*((i<0)?(int)Ny/Procy:(int)Ny/Procy + (int)Ny%Procy) ;
  
  
  
    disp_polymer[i]= i*((int)Nx/Procx)*((int)Ny/Procy);//recvcounts_polymer[i-1] + disp_polymer[i-1];    
   
   printf("Rank=%d::disp_poly[%d]=%d, count[%d]=%d\n",my2drank,i, disp_polymer[i],i,recvcounts_polymer[i]); 

   }
    //printf("Rank=%d,countx=%d, county=%d, total=%d\n", my2drank, count_localx1,count_localy1,count_localx1*count_localy1);



//--------------------------------------------------------------------------------------
  
     
      Nglobal_startx=3 + (int) Nx/Procx*coords[0]; 
      Nglobal_endx= Nglobal_startx + nlocalx-1;
   
      Nglobal_starty= 3+ (int) Ny/Procy*coords[1];
      Nglobal_endy=Nglobal_starty + nlocaly-1;
   
  
      printf("Rank=%d, (Nglobal_startx =%d, Nglobal_endx=%d)\n", my2drank, Nglobal_startx, Nglobal_endx);
	   
//------------------------------------------------------------------------------------

 

       PHI        = new double*[nlocalx+6];
       PHI_old    = new double*[nlocalx+6];
       gamma      = new double*[nlocalx+6];
       Laplacian2 = new double*[nlocalx+6];
       PP3        = new double[Nx*Ny];
       P2         = new double[Nx*Ny];

       //P2         = new double*[nlocalx+6];
       CP4        = new double*[nlocalx+6];
       PPP        = new double[Nx*Ny];
       P3         = new double*[nlocalx+6];
       PHI_local_result = new double[nlocalx*nlocaly];
       PPP_local        = new double[nlocalx*nlocaly];
       P2_local         = new double[nlocalx*nlocaly];
       POLYMER = new double [Nx*Ny]; //<----used for storing polymer points for coupling.
       matrix_lower = new double[3*nlocaly];
       matrix_upper = new double[3*nlocaly];
       matrix_left  = new double[3*(nlocalx+6)];
       matrix_right = new double[3*(nlocalx+6)];
      // std::shared_ptr<Node> nd( new Node);//malloc(sizeof(Node));

       for(int i=0;i<nlocalx+6; i++)
	  {
	   PHI[i]=new double[nlocaly+6];
	   PHI_old[i]=new double[nlocaly+6]; 
	   gamma[i]=new double[nlocaly+6];
	   Laplacian2[i]=new double [nlocaly+6];
	  // PP3[i]  =new double [nlocaly+6];
	   //P2[i]   = new double [nlocaly+6];
	   CP4[i]  = new double [nlocaly+6];
	   //PPP[i]  = new double [nlocaly+6];
	   P3[i]   = new double [nlocaly+6];


          }
      

       for(int i=0;i<nlocalx +6; i++)
	   {
            for(int j=0;j<nlocaly+6;j++)
	    {
               PHI[i][j]=0.0;
	       PHI_old[i][j]=0.0;
	       CP4[i][j]=0.0;
	    //   PP3[i][j]=0.0;
	       //P2[i][j]=0.0;
	       gamma[i][j]=0.0;
	       Laplacian2[i][j]=0.0;
	      // PPP[i][j]=0.0;
	       P3[i][j]=0.0;

	    }

	   }


//   XI.resize(nlocal_particles);
  // bd.resize(Nparticles);  
   bd = new BODY[Nparticles];
   Posx.resize(Nparticles);
   Posy.resize(Nparticles);
   Fx.resize(Nparticles);
   Fy.resize(Nparticles);
   Vx.resize(Nparticles);
   Vy.resize(Nparticles);
   //Recv_bodies= new BODY[nlocal_particles];
   //Send_bodies= new BODY[nlocal_particles];
   //down_nbr_list.resize(nlocal_particles);
   //right_nbr_list.resize(nlocal_particles);
   //up_nbr_list.resize(nlocal_particles);
   //left_nbr_list.resize(nlocal_particles);
   //Virtual_list.resize(nlocal_particles);  
   //Recv_bodies.resize(nlocal_particles);
   //send_bodies.resize(nlocal_particles);
   R12.resize(Nparticles);
   Dx.resize(Nparticles);
   Dy.resize(Nparticles);
   RI.resize(Nparticles);
   Dx1.resize(Nx);
   Dy1.resize(Nx);
   //PPP.resize(nlocalx);
   POS.resize(Nx);//<--- rethink about the size of this, it is bether to use a linked list
   VerletList_local.resize(Nparticles);
   VerletTable_local.resize(Nparticles);
   TOTX.resize(10*Nparticles);
   TOTY.resize(10*Nparticles);
   for(int i=0;i<Nx;i++)
      {
	POS[i].resize(Ny);
	Dx1[i].resize(Ny);
	Dy1[i].resize(Ny);
      }

   for(int i=0; i<Nparticles;i++)
      {
        Dx[i].resize(Nparticles);
        Dy[i].resize(Nparticles);
        R12[i].resize(Nparticles);
        RI[i].resize(Nparticles);
        VerletTable_local[i].resize(Nparticles);

        }
     


  up1x.resize(nlocaly+6);
  down1x.resize(nlocaly+6);


     printf("Base Construction made. \n");
  FreeCommunicator(new_comm); 

 }

BASE::~BASE()
     {
        for(int i=0;i<nlocalx+6;i++)
           {
             delete [] PHI[i];
	     delete [] PHI_old[i];
	     delete [] gamma[i];   
	     delete [] Laplacian2[i];
	     //delete [] PP3[i];
	   //  delete [] P2[i];
	     delete [] CP4[i];
	     //delete [] PPP[i];
	     delete [] P3[i];

	   }
	//for(int i=0;i<Nx;i++)
	  // {
             

	  // }
	  delete [] PHI;
	  delete [] PHI_old;     
	  delete [] Laplacian2;
	  delete [] gamma;	    
	  delete [] P3;
	  delete [] bd;
	 // delete [] Recv_bodies;
	 // delete [] Send_bodies;
	  delete [] recvcounts;
	  delete [] displacement;
	  delete [] recvcounts_polymer;
	  delete [] disp_polymer;
	  delete [] PP3;
	  delete [] P2;
	  delete [] CP4;
	  delete [] PPP;
	  delete [] PHI_local_result;
	  delete [] POLYMER;
	  delete [] PPP_local;
	  delete [] P2_local;

	  delete [] matrix_lower;
	  delete [] matrix_upper;
	  delete [] matrix_left;
	  delete [] matrix_right;
	 // delete  nd;
	  std::cout<<"Destructor in base class called"<<std::endl;
	//  std::cout<<"Destructor in base class called: colloids destroyed"<<std::endl; 

          //MPI_Comm_free(&new_comm);
	//  FreeCommunicator(new_comm);


     }
void BASE::InitializePolymer(MPI_Comm new_comm)
    {

      MPI_Comm_rank(new_comm, &my2drank);
      MPI_Comm_size(new_comm, &size);

     
      srand(time(NULL) + my2drank);
       for(int i=3; i<=nlocalx+2;i++)
            {
             for(int j=3;j<=nlocaly+2;j++)
	        {
	         double range =Random_max-Random_min;
	         double div =RAND_MAX / range;
	         PHI_old[i][j]=Random_min + (rand()/div);
		 PHI_local_result[(i-3)*nlocaly + j-3]= PHI_old[i][j];
                 PPP_local[(i-3)*nlocaly + j-3]=PPP_local[(i-3)*nlocaly + j-3] + B*pow((PHI_old[i][j]*P00b),2);
		 //printf("Rank=%d, PHI_old[%d][%d]=%lf\n", my2drank, i,j, PHI_old[i][j]);
	        }
	    }
    if(size>1)
     {
     MPI_Allgatherv(PHI_local_result, nlocalx*nlocaly,MPI_DOUBLE, POLYMER, recvcounts_polymer, disp_polymer,MPI_DOUBLE, new_comm );
     MPI_Allgatherv(PPP_local, nlocalx*nlocaly,MPI_DOUBLE, PP3, recvcounts_polymer, disp_polymer,MPI_DOUBLE, new_comm );
      }


     }

void BASE::setIndex(MPI_Comm new_comm)
{
   int coords[2], dims[2], periods[2];
   MPI_Comm_rank(new_comm, &my2drank);
   MPI_Comm_size(new_comm, &size);
   MPI_Cart_get(new_comm,2,dims,periods, coords );
   if(dims[1]==1)
   {
   for(int s=3; s<nlocaly+3; s++)
      {    
        up1x[s]=s+1;
        down1x[s]=s-1;
      }
    down1x[3]=nlocaly+2;
    up1x[nlocaly+2]=3;
   }
   else
   {
    for(int s=3; s<nlocaly+3; s++)
       {    
         up1x[s]=s+1;
         down1x[s]=s-1;
       }
   }
 }

void BASE::setLaplacianBase(MPI_Comm new_comm)
    {
      MPI_Comm_rank(new_comm, &my2drank);
      MPI_Comm_size(new_comm, &size);
         
         setIndex(new_comm);
         //cout<<"gamma(0,0)="<<gamma(0,0)<<"\t";                                                                       
         double AP=0.0, BP=0.0,  ATP=0.0;//gtemp=0.0;                                                           
         for(int i=3;i<=nlocalx+2;i++)
            {
            for(int j=3;j<=nlocaly+2;j++)
	       {
	        AP=(1.0/6.0)*(PHI_old(i+1,j) + PHI_old(i-1,j)	+ PHI_old(i,up1x[j]) + PHI_old(i,down1x[j]));
                BP=(1.0/12.0)*(PHI_old(i-1,up1x[j]) + PHI_old(i-1,down1x[j]) + PHI_old(i+1,up1x[j]) + PHI_old(i+1,down1x[j]));
                ATP = AP + BP;
	        gamma(i,j)=g(PHI_old(i,j))+ D*(ATP-PHI_old(i,j))-PHI_old(i,j);
               }
            }

         //polymer colloid coupling                                                                                
         double AP3(0.0), BP3(0.0);
         for(int i=3;i<=nlocalx+2;i++)
            {
             for(int j=3;j<=nlocaly+2;j++)
	       {
	        AP3=(1.0/6.0)*(P2(i+1,j)+P2(i-1,j) + P2(i,up1x[j]) + P2(i,down1x[j]));
	        BP3=(1.0/12.0)*(P2(i-1,down1x[j]) + P2(i-1,up1x[j]) + P2(i+1,down1x[j]) + P2(i+1,up1x[j]));
                CP4(i,j)=AP3 + BP3;
                //  printf("CP4(%d,%d)=%lf\n",i,j,CP4(i,j));                                                         
	       }
            }
   
 }

double BASE::g(double phi)
{
    double q(0.0);
    q=(1.0 + tau - A*pow((1.0-2.0*F),2))*phi-v*(1.0-2.0*F)*pow(phi,2)-u*pow(phi,3);
    return q;
}

void BASE::SetSecondLaplacian2(MPI_Comm new_comm)
{
      int coords[2], dims[2], periods[2];
      MPI_Comm_rank(new_comm, &my2drank);
      MPI_Comm_size(new_comm, &size);
      MPI_Cart_get(new_comm,2,dims,periods, coords );
              setIndex(new_comm);
          double AP=0.0, BP=0.0;
          for(int i=3;i<=nlocalx+2;i++)
             {
	      for(int j=3;j<=nlocaly+2;j++)
	         {
	           AP=(1.0/6.0)*(gamma(i+1,j)    + gamma(i-1,j) + gamma(i,up1x[j])    + gamma(i,down1x[j]));
	           BP=(1.0/12.0)*(gamma(i-1,up1x[j]) + gamma(i-1,down1x[j]) +gamma(i+1,up1x[j])  + gamma(i+1,down1x[j]));
	           Laplacian2(i,j) = AP + BP;
	         }
              }
}
void BASE::UpdateSolution(MPI_Comm new_comm)
{
      MPI_Comm_rank(new_comm, &my2drank);
      MPI_Comm_size(new_comm, &size);
      
          for(int i=3;i<nlocalx+3;i++)
	     {
	      for(int j=3;j<nlocaly+3;j++)
	         {
	          PHI_old[i][j]=PHI[i][j];
	          P3[i][j]=PHI_old[i][j];
		  PPP_local[(i-3)*nlocaly + j-3]=PPP_local[(i-3)*nlocaly + j-3] + B*pow((PHI_old[i][j]*P00b),2);
	         }
	     } 
     if(size>1)
      {
       MPI_Allgatherv(PHI_local_result, nlocalx*nlocaly,MPI_DOUBLE, POLYMER, recvcounts_polymer, disp_polymer,MPI_DOUBLE,new_comm);
       MPI_Allgatherv(PPP_local, nlocalx*nlocaly,MPI_DOUBLE, PP3, recvcounts_polymer, disp_polymer,MPI_DOUBLE,new_comm);
      }

  }

void BASE::FiniteDifferenceScheme(MPI_Comm new_comm)
{
      int coords[2];
      MPI_Comm_rank(new_comm, &my2drank);
      MPI_Comm_size(new_comm, &size);
      MPI_Cart_coords(new_comm,my2drank,2,coords);


      std::cout<<"in finite difference scheme parallel"<<std::endl;
     
	
         for(int i=3;i<nlocalx+3;i++)
            {
             for(int j=3;j<nlocaly+3;j++)
                {
		  //PHI[i][j]= PHI_old(i,j) - B*(PHI_old(i,j)-1.0 + 2*r)+ gamma(i,j)-Laplacian2(i,j);
	          int X=i+coords[0]*nlocalx; //<--global x index
		  int Y=j+coords[0]*nlocaly;
		  PHI[i][j]= PHI_old(i,j)+EMME*delta_t*(-B*(1.0-PPP(X,Y))*PHI_old(i,j)+ dxx*(gamma(i,j)-Laplacian2(i,j)) + dxx*CP4[i][j]);
	          PHI_local_result[(i-3)*nlocaly + j-3]= PHI[i][j];
	//	  printf("B=%lf\n",B);
                // printf("PHI_old[%d][%d]=%lf\n",i,j,PHI_old[i][j]); 
		}
             }

       
}

void BASE::ExchangeData(MPI_Comm new_comm, double **array)
{
     MPI_Status status;
     MPI_Comm_rank(new_comm,&my2drank);
     int coords[2], upper_nbr(0), down_nbr(0), tag, left_nbr(0), right_nbr(0) ;
     tag=201;
     int dims[2];
     int periods[2];
     MPI_Cart_coords(new_comm,my2drank,2,coords);
     MPI_Cart_get(new_comm,2,dims,periods, coords );
     //copy boundary points for exchange among mpi processes
    if((dims[0]>1)&&(dims[1]>1))
      {
       if(dims[0]>1)
        {
	   for(int i=3;i<6;i++)
	      {
                for(int j=3;j<=nlocaly+2; j++)
                   {
                    matrix_upper[(i-3)*nlocaly + j-3]=array[i][j];
		    //matrix_upper[j-1]=array[1][j];

                   }
	      } 
	   for(int i=nlocalx; i<nlocalx+3 ; i++)
	      {
               for(int j=3;j<=nlocaly+2; j++)
   	          {
                   matrix_lower[(i-nlocalx)*nlocaly + j-3]=array[i][j];
		   // matrix_lower[j-1]=array[nlocalx][j];
                  }
	      }
        MPI_Comm_rank(new_comm,&my2drank);
        MPI_Cart_coords(new_comm,my2drank,2,coords);
        MPI_Cart_shift( new_comm, 0, 1, &upper_nbr, &down_nbr ); //move along x direction
        MPI_Sendrecv_replace(matrix_upper, 3*nlocaly,  MPI_DOUBLE, upper_nbr, tag, down_nbr,tag, new_comm,&status);
        MPI_Sendrecv_replace(matrix_lower, 3*nlocaly,  MPI_DOUBLE, down_nbr, tag, upper_nbr,tag, new_comm,&status);
        //copy received ghost points into matrix, array
	for(int i=0;i<3;i++)
	   {
            for(int j=3; j<=nlocaly+2; j++)
               {
                array[i][j]=matrix_lower[(i)*nlocaly+j-3];
	       }
	   }
	for(int i=nlocalx+3; i<nlocalx+6; i++)
	   {
	     for(int j=3;j<nlocaly+3;j++)
		 {
                  array[i][j]=matrix_upper[(i-(nlocalx+3))*nlocaly + j-3];
		 }
	   }
     	   
      }
     if(dims[1]>0)
       {
         // for(int j=2;j<4;j++)
           //  {
	       for(int i=0;i<=nlocalx+5;i++)
                  {
	             matrix_left[i]=array[i][5];
                     matrix_left[(nlocalx+6)+i]=array[i][4];
		     matrix_left[2*(nlocalx+6)+i]=array[i][3];

		    }
             //   }
	  

        //for(int j=nlocaly;j<nlocaly+2;j++)
	  // {
             for(int i=0;i<=nlocalx+5;i++)
  	        {
		 matrix_right[ i]=array[i][nlocaly];
                 matrix_right[ (nlocalx+6)+i]=array[i][nlocaly+1];
		 matrix_right[2*(nlocalx+6) + i]=array[i][nlocaly+2];

		}
          // }
	 MPI_Comm_rank(new_comm,&my2drank);
	 MPI_Cart_coords(new_comm,my2drank,2,coords);
	 MPI_Cart_shift( new_comm, 1, 1, &left_nbr, &right_nbr );
	 MPI_Sendrecv_replace(matrix_right, 3*(nlocalx+6),  MPI_DOUBLE, right_nbr, tag, left_nbr,tag, new_comm,&status);
         MPI_Sendrecv_replace(matrix_left, 3*(nlocalx+6),  MPI_DOUBLE, left_nbr, tag, right_nbr,tag, new_comm,&status);
        //copy received ghost points into matrix
         
	    //for(int j=0;j<2;j++)
              // {
	        for(int i=0;i<=nlocalx+5; i++)
     	           {
         	     array[i][0]=matrix_right[i];
		     array[i][1]=matrix_right[(nlocalx+6)+i];
                     array[i][2]=matrix_right[2*(nlocalx+6)+i];
	           }
	      // }
	    //for(int j=nlocaly+3;j>(nlocaly+1);j--)
	      // {
	         for(int i=0;i<nlocalx+6;i++)
	            {
	             array[i][nlocaly+5]=matrix_left[i];
		     array[i][nlocaly+4]=matrix_left[(nlocalx+6)+i];
		     array[i][nlocaly+3]=matrix_left[2*(nlocalx+6)+i];

    	            }
	   
           //}
   }
 }
 else if((dims[0]>1)&&(dims[1]==1))
        {
         for(int j=3;j<=nlocaly+2; j++)
            {
              matrix_upper[ j-3]         =array[5][j];
	      matrix_upper[nlocaly + j-3]=array[4][j];
	      matrix_upper[2*nlocaly+ j-3]=array[3][j];
            }
	      //} 
	  // for(int i=nlocalx; i<nlocalx+2 ; i++)
	   //   {
               for(int j=3;j<=nlocaly+2; j++)
   	          {
                    matrix_lower[ j-3]         = array[nlocalx][j];
		    matrix_lower[nlocaly + j-3]= array[nlocalx+1][j];
                    matrix_lower[2*nlocaly+j-3]= array[nlocalx+2][j];
                  }
	     // }
        MPI_Comm_rank(new_comm,&my2drank);
        MPI_Cart_coords(new_comm,my2drank,2,coords);
        MPI_Cart_shift( new_comm, 0, 1, &upper_nbr, &down_nbr ); //move along x direction
        MPI_Sendrecv_replace(matrix_upper, 3*nlocaly,  MPI_DOUBLE, upper_nbr, tag, down_nbr,tag, new_comm,&status);
        MPI_Sendrecv_replace(matrix_lower, 3*nlocaly,  MPI_DOUBLE, down_nbr, tag, upper_nbr,tag, new_comm,&status);
        //copy received ghost points into matrix, array
	//for(int i=0;i<2;i++)
	  // {
            for(int j=3; j<=nlocaly+2; j++)
               {
                array[0][j]=matrix_lower[j-3];
		array[1][j]=matrix_lower[nlocaly+j-3];
		array[2][j]=matrix_lower[2*nlocaly+j-3];
	       }
	  // }
      //	for(int i=nlocalx+2; i<nlocalx+4; i++)
	//   {
	     for(int j=3;j<nlocaly+3;j++)
		 {
                   array[nlocalx+5][j] = matrix_upper[ j-3];
	           array[nlocalx+4][j] = matrix_upper[nlocaly + j-3];
		   array[nlocalx+3][j] = matrix_upper[2*nlocaly + j-3];
		   
		 }
	  // }
	
        }
 else if((dims[0]==1)&&(dims[1]>1))
       {
        std:: cout<<"Not yet implemented. Realigned the processors. That is, let Procy=1 and Procx=size"<<std::endl;
       }

    }



 

//-----------------END of POLYMER PART--------------------------------------------------------

//bool BASE:: IsZero(std::vector<particletype> *body)
//  {
  //   return (body->particletype[1]==0);//&&(bd[i].r[1]==0.0));
//  }

void BASE::ComputeForce(MPI_Comm new_comm)
{
   int coords[2],dims[2],periods[2];
   MPI_Cart_coords(new_comm,my2drank,2,coords);
   //MPI_Status status;
   //MPI_Request request;

   MPI_Comm_rank(new_comm, &my2drank);
   MPI_Comm_size(new_comm, &size);
   MPI_Cart_get(new_comm,2,dims,periods, coords);
   double G4(0.0), AA(0.0), BB(0.0), NEW1(0.0), NEW2(0.0), G2(0.0), G3(0.0), G3A(0.0),C(0.0);
   double DSUMX(0.0),DSUMY(0.0), DSUMXY(0.0);

for(int i=N_start;i<N_start+ nlocal_particles-1;i++)
     {
     for(int j=i+1;j<Nparticles;j++)
        {
               
	  double dx=bd[i].r[0]-bd[j].r[0];
	  double dy=bd[i].r[1]-bd[j].r[1];
	  //printf("dx=%lf, dy=%lf\n", dx,dy);
	 
	  double R2=sqrt(pow(dx,2)+pow(dy,2));
	  if(R2<2.0*R1big)
	    {
	     //std::cout<<"I am here"<<std::endl;
	     Dx[i][j]=dx;
	     Dy[i][j]=dy;
	     RI[i][j]=R2;
	     RI[j][i]=R2;
	     //C = U1/R12[i][j];
	     C = U1/(2.0*R1big);
	     //AA=RI[i][j]/R12[i][j];
	     AA=R2/(2.0*R1big);
	     G2 = pow(AA + beta,2.0);
	     G3 = 1.0 + ALPHA*(AA + beta);
	     BB = (AA-1.0)*(-ALPHA);
	     G3A = exp(BB);
	     NEW1=pow((1.0 + (ENNE/ALPHA)+ beta),2.0);
	     NEW2=(exp(-ENNE))*((1.0 + ALPHA*(beta + (ENNE/ALPHA)+1.0)))/(NEW1);
	     G4 = ((C*(G3A)*G3)/G2)-(C*NEW2);
	     bd[i].Force[0]=bd[i].Force[0] + G4*(dx/R2);//(Dx[i][j]/RI[i][j]);
	     bd[i].Force[1]=bd[i].Force[1] + G4*(dy/R2);//(Dy[i][j]/RI[i][j]);
	     //printf("Colloid Force (fx=%lf, fy=%lf) \n", bd[i].Force[0], bd[i].Force[1]);
	     bd[j].Force[0]=bd[j].Force[0] + G4*(Dx[j][i])/RI[i][j];
	     bd[j].Force[1]=bd[j].Force[1] + G4*(Dy[j][i])/RI[i][j];

	   }
	 }
   }
   //Compute random force
    double csi1(0.0), csi2(0.0);
    srand(time(NULL) + my2drank);
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
           bd[i].v[0]= (1.0/GAMM)*(bd[i].Force[0]+ bd[i].TOT[0])*delta_t + sqrt(2.0*temp1/GAMM)*csi1*sqrt(delta_t);
	   bd[i].v[1]= (1.0/GAMM)*(bd[i].Force[1]+ bd[i].TOT[1])*delta_t + sqrt(2.0*temp1/GAMM)*csi1*sqrt(delta_t);
	 // printf("Colloid Force (fx=%lf, fy=%lf) \n", Fx[i], Fy[i]);
	  }
	if(bd[i].index==1)
	  {
	   GAMM= GAMMb;
	   bd[i].v[0]= (1.0/GAMM)*(bd[i].Force[0]+ bd[i].TOT[0])*delta_t + sqrt(2.0*temp1/GAMM)*csi2*sqrt(delta_t);
	   bd[i].v[1]= (1.0/GAMM)*(bd[i].Force[1]+ bd[i].TOT[1])*delta_t + sqrt(2.0*temp1/GAMM)*csi2*sqrt(delta_t);
	   //printf("Colloid (TOTX =%lf, TOTY=%lf) \n", TOTX[i], TOTY[i]);
	   //printf("Colloid force (fx=%lf, fy=%lf) \n", bd[i].Force[0], bd[i].Force[1]);

           }
	   //position update
	  bd[i].r[0] =bd[i].r[0] + bd[i].v[0]; //<---bd[i].vx is already multiplied by t. 
	  bd[i].r[1] =bd[i].r[1] + bd[i].v[1];

          //bc condition
          
	  
	   bd[i].r[0]=bd[i].r[0]-Nx*( (int)(bd[i].r[0]/Nx + 1)-1);
	 
	   bd[i].r[1]=bd[i].r[1]-Ny*( (int)(bd[i].r[1]/Nx + 1)-1);
           //if(my2drank==0)	  
 	   //printf("Force_x[%d]=%lf, Force_y[%d]=%lf\n",i, bd[i].Force[0],i,bd[i].Force[1]);

	
        //computing mean displacement
        Posx[i]=Posx[i] + bd[i].v[0];
        Posy[i]=Posy[i] + bd[i].v[1];
        DSUMX  = DSUMX + pow(Posx[i],2.0);
        DSUMY  = DSUMY + pow(Posy[i],2.0);
        DSUMXY = DSUMXY + Posx[i]*Posy[i];
      //printf("Colloid positions (bd[100].x=%lf, bd[100].y=%lf) \n", bd[100].r[0], bd[100].r[1]);
    }
    // DSUMX=DSUMX/
}





void BASE::WriteToFile_MPI(MPI_Comm new_comm)
{
   MPI_Status status;
   MPI_File     fh;
   MPI_Datatype filetype;
   int dims[2],coords[2], periods[2], start_indices[2], localsizes[2];
   int globalsizes[2];
   globalsizes[0]=Nx;
   globalsizes[1]=Ny;
   localsizes[0]=(Nx/Procx); localsizes[1]=(Ny/Procy) ;
   int local_array_size =nlocalx*nlocaly;
   MPI_Cart_get (new_comm,2,dims,periods, coords );
   MPI_Comm_rank(new_comm, &my2drank);
   MPI_Comm_size(new_comm, &size);
   MPI_Cart_coords(new_comm,my2drank,2,coords);
   start_indices[0]=coords[0]*localsizes[0];
   start_indices[1]=coords[1]*localsizes[1];
   printf("rank=%d,r=%d,start_indices=(%d,%d)\n",my2drank,Nx%Procx,start_indices[0],start_indices[1]);
   MPI_Type_create_subarray(2, globalsizes, localsizes, start_indices,MPI_ORDER_C, MPI_DOUBLE, &filetype);
   MPI_Type_commit(&filetype);
   char outputfilename[]="datafile_mpi";
   char filememview[]="native";
   MPI_File_open(new_comm, outputfilename, MPI_MODE_CREATE | MPI_MODE_WRONLY,MPI_INFO_NULL, &fh);
   MPI_File_set_view(fh, 0, MPI_DOUBLE, filetype, filememview,MPI_INFO_NULL);
   MPI_File_write_all(fh, PHI_local_result,local_array_size , MPI_DOUBLE, &status);
   MPI_File_close(&fh);
   MPI_Type_free(&filetype);

}


void BASE::Read_input_parameters(int *integers, double *rationals)
     {
       //int value(0);
       FILE* file;
       char Data[MAX_LINE_LENGTH], *string;
	    if((file=fopen("ParameterFile.dat","r"))==NULL)
              {
	         printf("Error opening ParameterFile.dat\n");
	         return;
	      }
        
	string=fgets(Data,MAX_LINE_LENGTH,file);
          		
	int value=fscanf(file,"%d\n",&integers[0]);

	 string=fgets(Data,MAX_LINE_LENGTH,file);
         
	 value=fscanf(file,"%d\n",&integers[1]);
	 
         string=fgets(Data,MAX_LINE_LENGTH,file);
	 
	 value=fscanf(file,"%d\n",&integers[2]);
	 
	 string=fgets(Data,MAX_LINE_LENGTH,file);
         
	 value=fscanf(file,"%d\n",&integers[3]);
	 
	 string=fgets(Data,MAX_LINE_LENGTH,file);
	 
	 value=fscanf(file,"%d\n",&integers[4]);
	 
	string=fgets(Data,MAX_LINE_LENGTH,file);
         
	value=fscanf(file,"%d\n",&integers[5]);
	 
	string=fgets(Data,MAX_LINE_LENGTH,file);
	value=fscanf(file,"%d\n",&integers[6]);
        string=fgets(Data,MAX_LINE_LENGTH,file);
	value=fscanf(file,"%d\n",&integers[7]);
        string=fgets(Data,MAX_LINE_LENGTH,file);
 	value=fscanf(file,"%lf\n",&rationals[0]);
	string=fgets(Data,MAX_LINE_LENGTH,file);
        value=fscanf(file,"%lf\n",&rationals[1]);
	string=fgets(Data,MAX_LINE_LENGTH,file);
        value=fscanf(file,"%lf\n",&rationals[2]);
        string=fgets(Data,MAX_LINE_LENGTH,file);
        value=fscanf(file,"%lf\n",&rationals[3]);
        string=fgets(Data,MAX_LINE_LENGTH,file);
        value=fscanf(file,"%lf\n",&rationals[4]);
        string=fgets(Data,MAX_LINE_LENGTH,file);
        value=fscanf(file,"%lf\n",&rationals[5]);
        string=fgets(Data,MAX_LINE_LENGTH,file);
        value=fscanf(file,"%lf\n",&rationals[6]);
        string=fgets(Data,MAX_LINE_LENGTH,file);
        value=fscanf(file,"%lf\n",&rationals[7]);
        string=fgets(Data,MAX_LINE_LENGTH,file);
        value=fscanf(file,"%lf\n",&rationals[8]);
        string=fgets(Data,MAX_LINE_LENGTH,file);
	value=fscanf(file,"%lf\n",&rationals[9]);
	string=fgets(Data,MAX_LINE_LENGTH,file);
        value=fscanf(file,"%lf\n",&rationals[10]);
	string=fgets(Data,MAX_LINE_LENGTH,file);
	value=fscanf(file,"%lf\n",&rationals[11]);
	string=fgets(Data,MAX_LINE_LENGTH,file);
        value=fscanf(file,"%lf\n",&rationals[12]);
	string=fgets(Data,MAX_LINE_LENGTH,file);
        value=fscanf(file,"%lf\n",&rationals[13]);
	string=fgets(Data,MAX_LINE_LENGTH,file);
        value=fscanf(file,"%lf\n",&rationals[14]);
	string=fgets(Data,MAX_LINE_LENGTH,file);
        value=fscanf(file,"%lf\n",&rationals[15]);
	string=fgets(Data,MAX_LINE_LENGTH,file);
        value=fscanf(file,"%lf\n",&rationals[16]);

	string=fgets(Data,MAX_LINE_LENGTH,file);
        value=fscanf(file,"%lf\n",&rationals[17]);

	string=fgets(Data,MAX_LINE_LENGTH,file);
        value=fscanf(file,"%lf\n",&rationals[18]);
	string=fgets(Data,MAX_LINE_LENGTH,file);
        value=fscanf(file,"%lf\n",&rationals[19]);
	string=fgets(Data,MAX_LINE_LENGTH,file);
        value=fscanf(file,"%lf\n",&rationals[20]);
	string=fgets(Data,MAX_LINE_LENGTH,file);
        value=fscanf(file,"%lf\n",&rationals[21]);
	string=fgets(Data,MAX_LINE_LENGTH,file);
        value=fscanf(file,"%lf\n",&rationals[22]);
	string=fgets(Data,MAX_LINE_LENGTH,file);
        value=fscanf(file,"%lf\n",&rationals[23]);
	string=fgets(Data,MAX_LINE_LENGTH,file);
        value=fscanf(file,"%lf\n",&rationals[24]);
	string=fgets(Data,MAX_LINE_LENGTH,file);
        value=fscanf(file,"%lf\n",&rationals[25]);
//	string=fgets(Data,MAX_LINE_LENGTH,file);
//        value=fscanf(file,"%lf\n",&rationals[26]);
	
	
       // fgets(Data,MAX_LINE_LENGTH,file);
	//printf("value=%d,character=%s\n",value,string[0]);
	if((value!=0)||(string!=NULL))
	  {
           printf("failed to read in\n");
	  }
        fclose(file);
 
}



MPI_Comm BASE::CreatCartesianTopology()
      {
         MPI_Comm   new_comm;
         if(Procx*Procy!=size)
           {
             printf("Number of processors is not factorizable. That is Px*Py=%d, is not equal to the number of processors=%d\n ",Procx*Procy,size);
             printf("Hint: Check Parameter.dat file to change the number of cores\n");

		
	    exit(0);
	  }
        const int ROWS=0;const int COL=1;	     
        int periods[2];
        periods[0]=1; periods[1]=1;
        int coords[2];
        int dims[2];
        dims[ROWS]=Procx;
        dims[COL] =Procy;
        MPI_Comm_size( MPI_COMM_WORLD, &size );
        MPI_Comm_rank(MPI_COMM_WORLD,&rank);
        MPI_Cart_create( MPI_COMM_WORLD, 2,dims, periods, 0, &new_comm );
        MPI_Comm_rank(new_comm,&my2drank);
        MPI_Cart_coords(new_comm,my2drank,2,coords);
        MPI_Comm_size( new_comm, &size );

         return new_comm;

 }

void BASE::GenerateColloids(MPI_Comm new_comm)
{
   
    std::cout<<"----------------------- Generating Colloids--------------------"<<std::endl;
    MPI_Comm_rank(new_comm, &my2drank);
    MPI_Comm_size(new_comm, &size);
    //MPI_Status status;

    int coords[2], periods[2];
    int dims[2];
    MPI_Cart_coords(new_comm,my2drank,2,coords);
    MPI_Cart_get(new_comm,2,dims,periods, coords );

   
 
     srand(time(NULL) + my2drank);
     char fname[200];
     sprintf(fname, "solution%d.dat",my2drank );
     FILE *fp;
     fp=fopen(fname, "w");
       
         //double rangex=(double)(Nglobal_endx-Nglobal_startx);
	 //double rangey=(double)(Nglobal_endy-Nglobal_starty);
   
       	
           bd[N_start].r[0]= (((double)rand()/RAND_MAX))*(nlocalx)+ (double)coords[0]*(nlocalx); 
           bd[N_start].r[1]=  (((double)rand()/RAND_MAX))*nlocaly+ (double)coords[1]*(nlocaly);
		   //(double)coords[0]*nlocalx*(1.0/Nx) is the offset in x direction
            bd[N_start].v[0]=0.0;
            bd[N_start].v[1]=0.0;
            bd[N_start].index=1;
	    bd[N_start].TOT[0]=0.0;
            bd[N_start].TOT[1]=0.0;
	    bd[N_start].Force[0]=0.0;
	    bd[N_start].Force[1]=0.0;
	    
	    // printf("Colloid (rank=%d, dx=%lf, dy=%lf) \n",my2drank, bd[i].r[0], bd[i].r[1]);
        int i=0;       
       for(i=N_start+1; i<N_start+nlocal_particles; i++)
          {
label1:     bd[i].r[0]= (((double)rand()/RAND_MAX))*(nlocalx)+ (double)coords[0]*(nlocalx) ; 
	    bd[i].r[1]= (((double)rand()/RAND_MAX))*nlocaly+ (double)coords[1]*(nlocaly) ;
	    bd[i].v[0]=0.0;
            bd[i].v[1]=0.0;
            bd[i].index=1;
	    bd[i].TOT[0]=0.0;
	    bd[i].TOT[1]=0.0;
            bd[i].Force[0]=0.0;
	    bd[i].Force[1]=0.0;

           for(int j=N_start;j<i;j++)
	      {
               double RX=bd[i].r[0]-bd[j].r[0];
	       double RY=bd[i].r[1]-bd[j].r[1];
	       double R2=sqrt(pow(RX,2)+pow(RY,2));
	       if(R2<RR1big*2.0)//<---they should not occupy the same position
	         {
		   //std::cout<<"Going to label 1"<<std::endl;
		   goto label1;
		  }
	      }
	       fprintf(fp, "%d  %lf   %lf\n",i, bd[i].r[0], bd[i].r[1]);
	   }
       


 	 fclose(fp);
   
	 std::cout<<"---------------------End of Colloid Generation-----------------"<<std::endl;


}

/*
int BASE::OnBoundary( int ku,int i)
{
   int kd,kdd;
   kd = ku/2; // x(0)|y(1)|z(2) direction / 
   kdd = ku%2; // Lower(0)|higher(1) direction /
   if (kdd == 0)
      return bd[i].r[kd]< RCUT;
   else if (kdd==1)
      return al[kd]-RCUT < bd[i].r[kd];
   else
      return 0;
	  
}

int BASE::CheckBoundary(std::vector<BODY> body, int i)
{
  
 return 0;	

}*/

 void  BASE::Solver()
  {
    double t=0.0;
    int count =0;
    MPI_Comm new_comm;
    new_comm = CreatCartesianTopology();
    InitializePolymer(new_comm);
    GenerateColloids(new_comm);
    
    while(t<Max_time_iteration)
      {
        std::cout<<"in while loop solver"<<std::endl;
	
              
	 ExchangeData(new_comm,PHI_old);
         
	 Coupling(new_comm);
	 ParticleManager(new_comm);
	 if(t%10==0)
	 {
          ComputeForce( new_comm);
	 }
	 
	  //ExchangeData(new_comm,P2);
	  setLaplacianBase( new_comm);
	  ExchangeData(new_comm,gamma);
	  SetSecondLaplacian2(new_comm);
	 
	  //ExchangeData(new_comm,Laplacian2);
	  FiniteDifferenceScheme(new_comm);
	  UpdateSolution(new_comm);
	 
         
      t+=1.0;
      count++;
      printf("time=%lf, count=%d\n",t,count);
      }
   WriteToFile_MPI( new_comm);
 
 FreeCommunicator(new_comm);

 }




/*
void BASE::SetEffectiveRadius(MPI_Comm new_comm)
{
  // printf("IN EFFECTIVE RAD FUNCTION");
   MPI_Comm_rank(new_comm, &my2drank);
   MPI_Comm_size(new_comm, &size);
   

  for( auto i=N_start;i<N_start + nlocal_particles-1;i++)//(int)body->size()-1;i++)
      {
	for(size_t j=i+1; j<bd.size();j++) //Nparticles becos of verlet list
	   {
	   // if(((*body)[i].index==0)&&((*body)[j].index==0))
	     // {
	       //R12[i][j]=RR1small*2.0;

	      //}
	    if((bd[i].index==1)&&((bd[j].index==1)))
	      {
	       R12[i][j]=RR1big*2.0;
	      }
	   // if(((*body)[i].index==1)&&((*body)[j].index==0))
	    //{
	      //R12[i][j]=RR1big + RR1small;
	    //}
	    //if(((*body)[i].index==0)&&((*body)[j].index==1))
	    //{
	      //R12[i][j]=RR1big + RR1small;

	    //}
	 }
     }

}
*/
void BASE:: InsertList( ParticleList **root_list, ParticleList *i )
{

   i->next = *root_list;
   *root_list=i;
}


void BASE::ParticleManager(MPI_Comm new_comm)
{
   int coords[2];// send_offset, recv_offset;	
   //int dims[2];, pred, succ;
   //MPI_Status status;
   //MPI_Request send_request, recv_request;
   MPI_Comm_rank(new_comm, &my2drank);
   MPI_Comm_size(new_comm, &size);
   MPI_Cart_coords(new_comm,my2drank,2,coords);
   
 
   
    
   for (int i = 0; i < size; i++ )
     { 
        recvcounts[i] = (( i < (Nparticles%size) )? nlocal_particles + 1: nlocal_particles); 
	
        displacement[i] = ((int)(i*nlocal_particles) + ((i < Nparticles%size ) ? i : (Nparticles%size)));       
     }

   /*succ =(my2drank+1)%size;
   pred =(my2drank-1+size)%size;
   send_offset=N_start;//my2drank*nlocal_particles;
   recv_offset=((my2drank-1+size)%size)*recvcounts[pred];
   for(int i=0;i<size-1; i++)
      {
       MPI_Isend(&bd[send_offset], nlocal_particles, particletype, succ, 0, MPI_COMM_WORLD, &send_request);
       MPI_Irecv(&bd[displacement[pred]], recvcounts[pred], particletype, pred, 0, MPI_COMM_WORLD, &recv_request);
       send_offset=((my2drank-i-1+size)%size)*nlocal_particles;
       recv_offset=((my2drank-i-2+size)%size)*nlocal_particles;
       MPI_Wait(&send_request, &status);
       MPI_Wait(&recv_request, &status);
 
      }
   */
   
   
        //MPI_Comm_rank(new_comm,&my2drank);
       // MPI_Cart_coords(new_comm,my2drank,2,coords);
        //MPI_Cart_shift( new_comm, 1, 1, &left_nbr, &right_nbr ); //move along x direction

        //MPI_Isend(&bd[N_start], nlocal_particles,particletype, right_nbr, 0, new_comm, &send_request);
        //MPI_Irecv(&bd[displacement[left_nbr]],recvcounts[left_nbr],particletype, left_nbr, 0,new_comm, &recv_request);
       
       // MPI_Wait(&send_request, &status);
	//MPI_Wait(&recv_request, &status);
       //compute force here
	//MPI_Comm_rank(new_comm,&my2drank);
        //MPI_Cart_coords(new_comm,my2drank,2,coords);
        //MPI_Cart_shift( new_comm, 0, 1, &up_nbr, &down_nbr );
       
	//MPI_Isend(&bd[N_start], (nlocal_particles+recvcounts[left_nbr]),particletype, down_nbr, 0, new_comm, &send_request);
        //MPI_Irecv(&bd[displacement[left_nbr]],recvcounts[left_nbr],particletype, left_nbr, 0,new_comm, &recv_request);
     



   




   //for(int i=0;i<Virt)
   if (size>1)
   {
     MPI_Allgatherv(bd+N_start, nlocal_particles,particletype, bd, recvcounts, displacement,particletype, new_comm );
 

   }
   
   //Virtual_list.resize(1);
   
   //std::cout<<"Size On Each Core="<<<<std::endl;
//delete ParList;


 
}



void BASE::GenerateVerlet(MPI_Comm new_comm)
{
    //check this implementation. Check if the total number of particles per process is Nparticles. 
//Nparticles becos, after the call of MPI_AllgatherV, all processes should have Nparticles and generate only a local verlet list.
printf("IN VERLET FUNCTION");
 MPI_Comm_rank(new_comm, &my2drank);

 MPI_Comm_size(new_comm, &size);
    for(int i=N_start; i<N_start+nlocal_particles;i++)
       {
         VerletList_local[i]=0;
         for(int j=i+1;j<Nparticles;j++)
            {
             double dx=bd[i].r[0]-bd[j].r[0];
             double dy=bd[i].r[1]-bd[j].r[1];
             dx=dx-nlocalx*(int((dx/(nlocalx/2)+3)/2)-1);
             dy=dy-nlocaly*(int((dy/(nlocaly/2)+3)/2)-1);
             double R2=sqrt(pow(dx,2)+pow(dy,2));
             if(R2<R12[i][j])
               {
                VerletList_local[i]=VerletList_local[i] +1;
                VerletTable_local[i][VerletList_local[i]]=j;
	        Dx[i][j]=dx; //distance between body i and body j
	 //       Dx[j][i]=-dx;
	        Dy[i][j]=dy;
	   //     Dy[j][i]=-dy;
	        RI[i][j]=R2;
	     //   RI[j][i]=R2;
	      }
	    }
      }
	

}
void BASE:: FreeCommunicator(MPI_Comm new_comm)
{
         MPI_Comm_free(&new_comm);
		    

}
void BASE::Coupling(MPI_Comm new_comm)
{
   MPI_Comm_rank(new_comm, &my2drank);
   MPI_Comm_size(new_comm, &size);
   int coords[2], dims[2],periods[2];
   MPI_Cart_coords(new_comm,my2drank,2,coords);
   MPI_Cart_get(new_comm, 2,dims, periods,coords);
  // std::vector<auto> integration_index; //<--stores indices for integration
  char fname[200];
  sprintf(fname, "coupling%d.dat",my2drank );

  FILE *fp;
  fp=fopen(fname, "w");

   
   int N1(0), GT1(0), GT2(0), X, Y;
   double RX1(0.0),RY1(0.0), R2(0.0), B1(0.0),B2(0.0), D1(0.0), D2(0.0);
  for(int i=N_start;i<nlocal_particles + N_start;i++)
     {
       //int 
      if (bd[i].index==1 )//&&((bd[i].r[0]<=(double)Nglobal_endx)&&(bd[i].r[0]>=(double)Nglobal_startx)&&(bd[i].r[1]<=(double)Nglobal_endy)&&(bd[i].r[1]>=(double)Nglobal_starty))) //big particles
       {
	 /*here, we should consider only particles that are in the current process. We start by checking
	  weather the particles are within the MPI process.  */
	  
          P00=P00b;
          N1=(int)(RR1big+1.0);
	  bd[i].TOT[0]=0.0;
	  bd[i].TOT[1]=0.0;
           for(int j=-N1;j<N1;j++)
	      {
              for(int k=-N1;k<N1;k++)
  		  {
                    GT1=(int)bd[i].r[0] + j;
		   
		    GT2=(int)bd[i].r[1] + k;
		   // GT2=(int)GT2;
		   // X=(int)GT1;
		    if(GT1>Nx-1)
		      {
                        GT1=GT1-(Nx-1)*((int)(GT1/Nx +1)-1);
		      }
		    if(GT1<0)
		      {
                        GT1= GT1 + Nx-1;
		      }
		    X = (int)GT1;

                    if(GT2>Ny-1)
		      {
                        GT2=GT2-(Ny-1)*((int)(GT2/Ny +1)-1);
		      }
		    if(GT2<0)
		      {
                        GT2= GT2 + Ny-1;
		      }
		    Y=(int)GT2;
		    
		  
                    RX1=bd[i].r[0]-(double)GT1;
		    RX1=RX1-(double)(Nx*((int)((RX1/(Nx/2)+3)/2)-1));
         	    RY1=bd[i].r[1]-(double)GT2;
		    RY1=RY1-(double)(Ny*((int)((RY1/(Ny/2)+3)/2)-1));

		  
		    R2=sqrt(pow(RX1,2)+pow(RY1,2));
		    if(R2<RR1big)
	              {
	               Dx1[X][Y]=RX1;
	               Dy1[X][Y]=RY1;
	               B1 = 1.0-pow(R2/RR1big, ALP1);
	               B2 =exp(1.0-(1.0/B1));
	               D1 =ALP1/pow(RR1big,ALP1);
	               D2 =-(1.0)*(D1*pow(R2,ALP1-2.0)/(pow(B1,2)))*B2;//<--check this
		       //printf("X=%d, Y=%d\n", X,Y);
		       POS[X][Y]=i;
		       //PP3(X,Y)=PP3(X,Y) + B2*(POLYMER(X,Y)-P00);

	               
		       //printf("Rank=%d, PP3[%d][%d]=%lf\n",my2drank, X, Y, PP3[X][Y]);
		       if(R2<R1big)
			 {
                           PPP(X,Y)=1.0;
			 }
	              bd[i].TOT[0]+= -sigma*RX1*D2*pow((POLYMER(X,Y)-P00),2);
	              bd[i].TOT[1]+= -sigma*RY1*D2*pow((POLYMER(X,Y)-P00),2);
		       fprintf(fp, "%d  %lf   %lf\n",i, bd[i].TOT[0], bd[i].TOT[1]);

		      // printf("Rank=%d,  (X =%d, Y=%d,P3[%d][%d]=%lf) \n",my2drank,X, Y,X,Y, POLYMER(X,Y));
		       //printf("Colloid (X =%d, Y=%d,TOTX=%lf,TOTY=%lf) \n",X, Y, bd[i].TOT[0],bd[i].TOT[1]);

	              }
       		   }
                }
              }
           } //end of particle loop
       fclose(fp);
   //computation equation 29
 // if(size>1)
   // {
     // MPI_Allgatherv(PHI_local_result, nlocalx*nlocaly,MPI_DOUBLE, POLYMER, recvcounts_polymer, disp_polymer,     MPI_DOUBLE, new_comm );
   // }

     
   for(int i=0;i<nlocalx;i++)
      {
       for(int j=0;j<nlocaly;j++)
	  {
           if(POS[i][j]!=0)
            {
             int z= POS[i][j];
	      if(bd[z].index==0)
		 {
                  P00=P00s;
		 }
	      else if(bd[z].index==1)
		 {
                   P00=P00b;
		 }
	      else
		 {
	          continue;
		 }
	    }
	   int glob_index_x = i + coords[0]*nlocalx;
	   int glob_index_y = j + coords[1]*nlocaly;
           P2_local[i*nlocalx+j]=sigma*PP3(glob_index_x,glob_index_y);
	   printf("PP3[%d][%d]=%lf\n", glob_index_x,glob_index_y, PP3(glob_index_x,glob_index_y));
	  }
      }
   if(size>1)
   {
     MPI_Allgatherv(P2_local, nlocalx*nlocaly,MPI_DOUBLE, P2, recvcounts_polymer, disp_polymer,MPI_DOUBLE,new_comm);
   }

   //calculate laplacian for use with polymer
   
   double AP3(0.0), BP3(0.0);

   for(int i=3;i<=nlocalx+2;i++)
      {
       for(int j=3;j<=nlocaly+2;j++)
	  {
            int glob_index_x = i-3 + coords[0]*nlocalx;
	    int glob_index_y = j-3 + coords[1]*nlocaly;

            AP3=(1.0/6.0)*(P2(glob_index_x+1,glob_index_y )+P2(glob_index_x -1,glob_index_y ) + P2(glob_index_x ,glob_index_y +1) + P2(glob_index_x ,glob_index_y-1));
	    BP3=(1.0/12.0)*(P2(glob_index_x-1,glob_index_y-1) + P2(glob_index_x-1,glob_index_y+1) + P2(glob_index_x+1, glob_index_y-1) + P2(glob_index_x+1,glob_index_y+1));
	    CP4(i,j)=AP3 + BP3;
	    //printf("CP4[%d][%d]=%lf\n", i,j, CP4(i,j));
	  }

      }
     
}

