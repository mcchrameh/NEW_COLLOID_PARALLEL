  int n_to_move_right=0, n_to_move_left=0, n_to_move_up=0, n_to_move_down=0;
  std::vector<BODY> newparticles;
  std::vector<MPI_Aint> elmoffset_down;
  elmoffset_down.resize(Nparticles);
  std::vector<int> elmsize_down;
  elmsize_down.resize(Nparticles);
  std::vector<MPI_Aint> elmoffset_right;
  std::vector<int> elmsize_right;
  elmoffset_right.resize(Nparticles);
  elmsize_right.resize(Nparticles);
 
  std::vector<MPI_Aint> elmoffset_left;
  std::vector<int> elmsize_left;
  elmoffset_left.resize(Nparticles);
  elmsize_left.resize(Nparticles);
    
  std::vector<MPI_Aint> elmoffset_up;
  std::vector<int> elmsize_up;
  elmoffset_up.resize(Nparticles);
 
  elmsize_up.resize(Nparticles);
  
  MPI_Datatype sendtype_down, sendtype_up, sendtype_left, sendtype_right;
  for(int i=N_start; i<N_start + nlocal_particles;i++)
       {
         if(bd[i].r[0]>(double)(Nglobal_endx))//if particles crosses the boundary
            {
              MPI_Address(&bd[N_start], &elmoffset_down[n_to_move_down]);
              elmsize_down[n_to_move_down]=1;
              n_to_move_down++;
              bd.erase(bd.begin()+N_start+i);
              printf("down_move=%d\n",n_to_move_down); 
           }
        else if(bd[i].r[1]>(double)(Nglobal_endy))
           {
             MPI_Address(&bd[N_start],&elmoffset_right[n_to_move_right]);
             elmsize_right[n_to_move_right]=1;
             n_to_move_right++;
             bd.erase(bd.begin()+N_start+i);
             printf("right_move=%d\n",n_to_move_right);
           }
       else if(bd[i].r[0]<(double)Nglobal_startx)
         {
            MPI_Address(&bd[N_start], &elmoffset_up[n_to_move_up]);
            elmsize_up[n_to_move_up]=1;
            bd.erase(bd.begin()+N_start+i);
            n_to_move_up++;
           }
       else if(bd[i].r[1]<(double)Nglobal_starty)
       {
         MPI_Address(&bd[N_start], &elmoffset_left[n_to_move_left]);
         elmsize_left[n_to_move_left]=1;
         bd.erase(bd.begin()+N_start+i);
         n_to_move_left++;
       }
     }
    
    if(n_to_move_down<Nparticles)
     {
        elmoffset_down.resize(n_to_move_down);
        elmsize_down.resize(n_to_move_down);
     }
   if(n_to_move_up<Nparticles)
      {
       elmoffset_up.resize(n_to_move_down);
       elmsize_up.resize(n_to_move_down);
      }
   if(n_to_move_left<Nparticles)
      {
        elmoffset_left.resize(n_to_move_down);
        elmsize_left.resize(n_to_move_down);
      }
  if(n_to_move_right<Nparticles)
      {
        elmoffset_right.resize(n_to_move_down);
        elmsize_right.resize(n_to_move_down);
      }
  MPI_Type_hindexed(n_to_move_down,&elmsize_down[0], &elmoffset_down[0], particletype, &sendtype_down );
  MPI_Type_hindexed(n_to_move_up, &elmsize_up[0], &elmoffset_up[0], particletype, &sendtype_up );
  MPI_Type_hindexed(n_to_move_left,&elmsize_left[0], &elmoffset_left[0], particletype, &sendtype_left );
  MPI_Type_hindexed(n_to_move_right,&elmsize_right[0], &elmoffset_right[0], particletype, &sendtype_right );
  MPI_Type_commit(&sendtype_down);
  MPI_Type_commit(&sendtype_up);
  MPI_Type_commit(&sendtype_left);
  MPI_Type_commit(&sendtype_right);
  int number=0;
  MPI_Cart_coords(new_comm,my2drank,2,coords);
  MPI_Cart_shift( new_comm, 0, 1, &upper_nbr, &down_nbr );
  MPI_Send(MPI_BOTTOM,1, sendtype_up, upper_nbr, tag, new_comm);
  MPI_Probe(down_nbr,tag, new_comm, &status);
  MPI_Get_count(&status,particletype, &number);
  MPI_Type_extent(particletype,&extent);
  newparticles.resize(number*extent);
  MPI_Recv(&newparticles[0], number, particletype,down_nbr, tag, new_comm, &status);
  MPI_Send(MPI_BOTTOM,1, sendtype_down, down_nbr, tag, new_comm);
  printf("number=%d\n",number);
  number =0;  //<----reset number
  MPI_Probe(upper_nbr,tag, new_comm, &status);
  MPI_Get_count(&status,particletype, &number);
  MPI_Type_extent(particletype,&extent);
  newparticles.resize(number*extent);
  MPI_Recv(&newparticles[0], number, particletype,upper_nbr, tag, new_comm, &status);
  //MPI_Irecv(&bd[N_start], n_to_move_up, particletype,upper_nbr, tag, new_comm, &request);
  MPI_Cancel(&request);
  MPI_Cart_coords(new_comm,my2drank,2,coords);
  MPI_Cart_shift( new_comm, 1, 1, &left_nbr, &right_nbr);
  //MPI_Send(&bd[N_start],1, sendtype_left, left_nbr, tag, new_comm);
 //MPI_Irecv(&bd[N_start], n_to_move_right, particletype,right_nbr, tag, new_comm, &request);
 MPI_Type_commit(&sendtype_down);
 MPI_Type_commit(&sendtype_up);
//MPI_Send(&bd[N_start],1, sendtype_right, right_nbr, tag, new_comm);
   //MPI_Irecv(&bd[N_start], n_to_move_left, particletype,left_nbr, tag, new_comm, &request);
 MPI_Type_free(&sendtype_left);
 MPI_Type_free(&sendtype_right);
 MPI_Type_free(&sendtype_up);
 MPI_Type_free(&sendtype_down);


