
   
   
void BASE::GenerateColloids(MPI_Comm new_comm)
 {
  
   std::cout<<"----------------------- Generating Colloids--------------------"<<std::endl;
    MPI_Comm_rank(new_comm, &my2drank);
    MPI_Comm_size(new_comm, &size);
 
     int coords[2], periods[2];
     int dims[2];
     MPI_Cart_coords(new_comm,my2drank,2,coords);
     MPI_Cart_get(new_comm,2,dims,periods, coords );
 
 
 
      srand(time(NULL) + my2drank);
      char fname[200];
      sprintf(fname, "solution%d.dat",my2drank );
     FILE *fp;
     fp=fopen(fname, "w");
     int i=0;
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
                    // printf("Colloid (rank=%d, dx=%lf, dy=%lf) \n",my2drank, bd[i].r[0], bd[i].r[1]);
                  fprintf(fp, "%d  %lf   %lf\n",N_start, bd[N_start].r[0], bd[N_start].r[1]);
        for(i=N_start+1; i<nlocal_particles+N_start; i++)
                  {
 label1:            bd[i].r[0]= (((double)rand()/RAND_MAX))*nlocalx + (double)coords[0]*(nlocalx);
                    bd[i].r[1]= (((double)rand()/RAND_MAX))*nlocaly +(double)coords[1]*(nlocaly);
 
                    bd[i].v[0]=0.0;
                    bd[i].v[1]=0.0;
                    bd[i].index=1;
                    bd[i].TOT[0]=0.0;
                    bd[i].TOT[1]=0.0;
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
       fclose(fp);
 
          std::cout<<"---------------------End of Colloid Generation-----------------"<<std::endl;
 
 
 }
		 
                                                                         
                                                                                             
