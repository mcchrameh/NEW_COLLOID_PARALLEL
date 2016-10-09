




/*
   void BASE::SendBodies( MPI_Comm new_comm)
   {

      MPI_Status status;
         MPI_Request request;	   
	    // ku represents the four corners/neighbors of a process /
	       int ku,ku_0,ku_1,nsd;
	          MPI_Comm_rank(new_comm, &my2drank);
		     MPI_Comm_size(new_comm, &size);
		        int N_recv=0, tag=1;
			   for(int kd=0;kd<2;kd++) //loop over x and y direction
			         {
				        for (int bc=0;bc<2;bc++)  //bc=0 for lower(left) boundary, bc=1 for upper(right) boundary 
						    {
						                  	    
						    	      lsb[2*kd+bc][0]=0; //reset the number of to be copied atoms on all 4 boundaries
							                  }
									         // scan all atoms both resident and copies to identify boundary atoms /
										        for(int i=0;i<nlocal_particles+N_recv;i++)
												   {
												                //MPI_Cart_shift(new_comm, kd, 1,ku_0, ku_1);
														             for(int bc=0;bc<2;bc++)
															     		{
																		          ku=2*kd + bc;
																			  	         if (OnBoundary(ku,i))
																					 		    {
																							                          lsb[ku][++(lsb[ku][0])]=i;
																										  		     // sendBodies.push_back(bd[i]);
																												     		    }
																														    	          }
																																            
																																               }
																																	             }
																																		              printf("Number of Boundary atoms,left=%d\n",lsb[0][0] );
																																			              //left boundary buffer
																																				               
																																				      	 send_bodies_left.resize(lsb[0][0]);
																																					 	 for(int i=0;i<lsb[0][0];i++)
																																						             {
																																							     	     send_bodies_left.push_back(bd[i]);
																																								     	    }
																																									    	 //right boundary buffer
																																										 	 send_bodies_right.resize(lsb[2][0]);
																																											 	 for(int i=0;i<lsb[2][0];i++)
																																												 	    {       
																																													    	     send_bodies_right.push_back(bd[i]);
																																														                 }
																																																 	 //upper boundary buffer
																																																	 	 send_bodies_upper.resize(lsb[1][0]);
																																																		 	 for(int i=0;i<lsb[1][0];i++)
																																																			 	    {       
																																																				                   send_bodies_upper.push_back(bd[i]);
																																																						               }
																																																							       	 //lower boundary buffer
																																																								 	 send_bodies_lower.resize(lsb[3][0]);
																																																									 	 for(int i=0;i<lsb[3][0];i++)
																																																										 	    {       
																																																											    	      send_bodies_lower.push_back(bd[i]);
																																																												                  }

																																																														             //get neighbors in x and y direction to send and recieve
																																																															     	    MPI_Comm_rank(new_comm,&my2drank);
																																																																    	    MPI_Cart_shift(new_comm,0, 1,&ku_0, &ku_1);
																																																																	    	   // MPI_Sendrecv_replace(&send_bodies_upper, lsb[1][0]*sizeof(BODY)/sizeof(double),  MPI_DOUBLE, ku_1, tag,
																																																																		   	   //                       ku_0,tag, new_comm,&status);
																																																																			   	    MPI_Isend(&send_bodies_upper,lsb[1][0]*sizeof(BODY)/sizeof(double),MPI_DOUBLE,ku_1,tag,new_comm,&request);
																																																																				    	    MPI_Irecv(&send_bodies_upper,lsb[1][0]*sizeof(BODY)/sizeof(double),MPI_DOUBLE,ku_0,tag,new_comm,&request);
																																																																					                int new_size=bd.size() + lsb[1][0];
																																																																								    bd.resize(new_size);
																																																																								    	    bd.insert(std::end(bd),std::begin(send_bodies_upper),std::end(send_bodies_upper) );

	    //MPI_Sendrecv_replace(&send_bodies_lower, lsb[3][0]*sizeof(BODY)/sizeof(double),  MPI_DOUBLE, ku_0, tag,
	   //	                    ku_1,tag, new_comm,&status);

	    MPI_Isend(&send_bodies_lower,lsb[3][0]*sizeof(BODY)/sizeof(double),MPI_DOUBLE,ku_0,tag,new_comm,&request);
	                MPI_Irecv(&send_bodies_lower,lsb[3][0]*sizeof(BODY)/sizeof(double),MPI_DOUBLE,ku_1,tag,new_comm,&request);
			            int new_size1=bd.size() + lsb[3][0];
				    	    
				    	    bd.resize(new_size1);
					    	    bd.insert(std::end(bd),std::begin(send_bodies_lower),std::end(send_bodies_lower) );

						    	    MPI_Comm_rank(new_comm,&my2drank);
							    	    MPI_Cart_shift(new_comm,1, 1,&ku_0, &ku_1);
								    	    //MPI_Sendrecv_replace(&send_bodies_left, lsb[0][0]*sizeof(BODY)/sizeof(double),  MPI_DOUBLE, ku_0, tag,
								    	    //                         ku_1,tag, new_comm,&status);
								    	    MPI_Isend(&send_bodies_left,lsb[0][0]*sizeof(BODY)/sizeof(double),MPI_DOUBLE,ku_0,tag,new_comm,&request);
									    	    MPI_Irecv(&send_bodies_left,lsb[0][0]*sizeof(BODY)/sizeof(double),MPI_DOUBLE,ku_1,tag,new_comm,&request);
										                int new_size2=bd.size() + lsb[0][0];   
												            bd.resize(new_size2);
													    	    bd.insert(std::end(bd),std::begin(send_bodies_left),std::end(send_bodies_left) );
														    	                
														    	   // MPI_Sendrecv_replace(&send_bodies_right, lsb[2][0]*sizeof(BODY)/sizeof(double),  MPI_DOUBLE, ku_1, tag,
														    	    //                        ku_0,tag, new_comm,&status);
														    	    MPI_Isend(&send_bodies_right,lsb[2][0]*sizeof(BODY)/sizeof(double),MPI_DOUBLE,ku_1,tag,new_comm,&request);
															    	    MPI_Irecv(&send_bodies_right,lsb[2][0]*sizeof(BODY)/sizeof(double),MPI_DOUBLE,ku_0,tag,new_comm,&request);

																    	    int new_size3=bd.size() + lsb[2][0];   
																	                bd.resize(new_size3);
																				    bd.insert(std::end(bd),std::begin(send_bodies_right),std::end(send_bodies_right) );

																				          
																				    }

*/
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
/*
   void BASE::DistributeColloids(MPI_Comm new_comm)
   {
      //MPI_Status status;
        
        MPI_Comm_rank(new_comm, &my2drank);
	  MPI_Comm_size(new_comm, &size);
	    std::vector<BODY>:: iterator it=bd.begin();
	     //use cyclic distribution
	       if((size>1)&&(it!=bd.end()))
	          {
		      // std::vector<BODY>::iterator it=bd.begin();
		           MPI_Allgatherv(bd + *it, nlocal_particles*sizeof(BODY)/sizeof(double),MPI_DOUBLE, bd, recvcounts, displacement,MPI_DOUBLE, new_comm );
			      }
			      }
			      */


