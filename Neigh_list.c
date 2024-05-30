#include<stdio.h>
#include<stdlib.h>
#include<math.h>

void Neigh_list(int nx,int ny,int nz,int NDX,int NDY,int NDZ,int Totaldomains,int Totalatoms,int Neigh[Totaldomains][27],int BC)
{
int Ix,Iy,Iz,N_Ix,N_Iy,N_Iz,E,N,p,q,r;
int Neighbor;
E=1;
if(BC==1)                      //[1,1,1] periodic boundary condition
{
for(Iz=1;Iz<=NDZ;Iz++)
{
 for(Iy=1;Iy<=NDY;Iy++)
  {
   for(Ix=1;Ix<=NDX;Ix++)
   {
    E=Iy+NDY*(Ix-1+NDX*(Iz-1));               //Cell_Id of the cell whose neighbors are to be found
     N=0;
    for(p=Iz-1;p<=Iz+1;p++)
     {
     for(q=Iy-1;q<=Iy+1;q++)
      {
      for(r=Ix-1;r<=Ix+1;r++)
       {
       N_Iz=p;                             //N_Ix,N_Iy,N_Iz are the Ix,Iy,Iz values of neighbor cell 
       N_Iy=q;
       N_Ix=r;
      
       if(N_Ix==0)
       N_Ix=NDX;
       else if(N_Ix>NDX)
       N_Ix=1;
       if(N_Iy==0)
       N_Iy=NDY;
       else if(N_Iy>NDY)
       N_Iy=1;
       if(N_Iz==0)
       N_Iz=NDZ;
       else if(N_Iz>NDZ)
       N_Iz=1;
      
       Neighbor= N_Iy+NDY*(N_Ix-1+NDX*(N_Iz-1));  //Id of neighbor cell
       Neigh[E][N]=Neighbor;                      //array containing neighbors of 'E'th cell.
       N++;
       }
      }
     }
    }
   }
  }
}

if(BC==2)                        //[1,1,0] boundary condition
{for(Iz=1;Iz<=NDZ;Iz++)
{
 for(Iy=1;Iy<=NDY;Iy++)
  {
   for(Ix=1;Ix<=NDX;Ix++)
   {
    E=Iy+NDY*(Ix-1+NDX*(Iz-1));
     N=0;
    for(p=Iz-1;p<=Iz+1;p++)
     {
     for(q=Iy-1;q<=Iy+1;q++)
      {
      for(r=Ix-1;r<=Ix+1;r++)
       {
       N_Iz=p;
       N_Iy=q;
       N_Ix=r;
      
       if(N_Ix==0)
       N_Ix=NDX;
       else if(N_Ix>NDX)
       N_Ix=1;
       if(N_Iy==0)
       N_Iy=NDY;
       else if(N_Iy>NDY)
       N_Iy=1;
       if(N_Iz==0)
       {N_Iz=NDZ;
       break;
       }
       else if(N_Iz>NDZ)
       {N_Iz=1;
       break;}
      else if(N_Iz!=0 || N_Iz<NDZ)
      {
       Neighbor= N_Iy+NDY*(N_Ix-1+NDX*(N_Iz-1));
       Neigh[E][N]=Neighbor;
       N++;
      }
       }
      }
     }
    }
   }
  }
}

if(BC==3)                           //[0,1,1] boundary condition
{for(Iz=1;Iz<=NDZ;Iz++)
{
 for(Iy=1;Iy<=NDY;Iy++)
  {
   for(Ix=1;Ix<=NDX;Ix++)
   {
    E=Iy+NDY*(Ix-1+NDX*(Iz-1));
     N=0;
    for(p=Iz-1;p<=Iz+1;p++)
     {
     for(q=Iy-1;q<=Iy+1;q++)
      {
      for(r=Ix-1;r<=Ix+1;r++)
       {
       N_Iz=p;
       N_Iy=q;
       N_Ix=r;
      
       if(N_Ix==0)
       {N_Ix=NDX;
        break;}
       else if(N_Ix>NDX)
       {N_Ix=1;
        break;}
       if(N_Iy==0)
       N_Iy=NDY;
       else if(N_Iy>NDY)
       N_Iy=1;
       if(N_Iz==0)
       {N_Iz=NDZ;
       break;
       }
       else if(N_Iz>NDZ)
       {N_Iz=1;
       break;}
      else if(N_Ix!=0 || N_Ix<NDX)
      {
       Neighbor= N_Iy+NDY*(N_Ix-1+NDX*(N_Iz-1));
       Neigh[E][N]=Neighbor;
       N++;
      }
       }
      }
     }
    }
   }
  }
}

if(BC==4)                         //[1,0,1] boundary condition
{for(Iz=1;Iz<=NDZ;Iz++)
{
 for(Iy=1;Iy<=NDY;Iy++)
  {
   for(Ix=1;Ix<=NDX;Ix++)
   {
    E=Iy+NDY*(Ix-1+NDX*(Iz-1));
     N=0;
    for(p=Iz-1;p<=Iz+1;p++)
     {
     for(q=Iy-1;q<=Iy+1;q++)
      {
      for(r=Ix-1;r<=Ix+1;r++)
       {
       N_Iz=p;
       N_Iy=q;
       N_Ix=r;
      
       if(N_Ix==0)
       {N_Ix=NDX;
        }
       else if(N_Ix>NDX)
       {N_Ix=1;
        }
       if(N_Iy==0)
       {N_Iy=NDY;
       break;}
       else if(N_Iy>NDY)
       {N_Iy=1;
        break;}
       if(N_Iz==0)
       {N_Iz=NDZ;
              }
       else if(N_Iz>NDZ)
       {N_Iz=1;
       }
      else if(N_Iy!=0 || N_Ix<NDY)
      {
       Neighbor= N_Iy+NDY*(N_Ix-1+NDX*(N_Iz-1));
       Neigh[E][N]=Neighbor;
       N++;
      }
       }
      }
     }
    }
   }
  }
}

if(BC==5)                 //[1,0,0] boundary condition
{for(Iz=1;Iz<=NDZ;Iz++)
{
 for(Iy=1;Iy<=NDY;Iy++)
  {
   for(Ix=1;Ix<=NDX;Ix++)
   {
    E=Iy+NDY*(Ix-1+NDX*(Iz-1));
     N=0;
    for(p=Iz-1;p<=Iz+1;p++)
     {
     for(q=Iy-1;q<=Iy+1;q++)
      {
      for(r=Ix-1;r<=Ix+1;r++)
       {
       N_Iz=p;
       N_Iy=q;
       N_Ix=r;
      
       if(N_Ix==0)
       {N_Ix=NDX;
        }
       else if(N_Ix>NDX)
       {N_Ix=1;
        }
       if(N_Iy==0)
       {N_Iy=NDY;
       break;}
       else if(N_Iy>NDY)
       {N_Iy=1;
        break;}
       if(N_Iz==0)
       {N_Iz=NDZ;
       break;
       }
       else if(N_Iz>NDZ)
       {N_Iz=1;
       break;}
      else if((N_Iy!=0 || N_Iy<NDY) && (N_Iz!=0 || N_Iz<NDZ))
      {
       Neighbor= N_Iy+NDY*(N_Ix-1+NDX*(N_Iz-1));
       Neigh[E][N]=Neighbor;
       N++;
      }
       }
      }
     }
    }
   }
  }
}

if(BC==6)                                 //[0,1,0] boundary condition
{for(Iz=1;Iz<=NDZ;Iz++)
{
 for(Iy=1;Iy<=NDY;Iy++)
  {
   for(Ix=1;Ix<=NDX;Ix++)
   {
    E=Iy+NDY*(Ix-1+NDX*(Iz-1));
     N=0;
    for(p=Iz-1;p<=Iz+1;p++)
     {
     for(q=Iy-1;q<=Iy+1;q++)
      {
      for(r=Ix-1;r<=Ix+1;r++)
       {
       N_Iz=p;
       N_Iy=q;
       N_Ix=r;
      
       if(N_Ix==0)
       {N_Ix=NDX;
        break;}
       else if(N_Ix>NDX)
       {N_Ix=1;
        break;}
       if(N_Iy==0)
       {N_Iy=NDY;
       }
       else if(N_Iy>NDY)
       {N_Iy=1;
       }
       if(N_Iz==0)
       {N_Iz=NDZ;
       break;
       }
       else if(N_Iz>NDZ)
       {N_Iz=1;
       break;}
      else if((N_Ix!=0 || N_Ix<NDX) && (N_Iz!=0 || N_Iz<NDZ))
      {
       Neighbor= N_Iy+NDY*(N_Ix-1+NDX*(N_Iz-1));
       Neigh[E][N]=Neighbor;
       N++;
      }
       }
      }
     }
    }
   }
  }
}

if(BC==7)                                 //[0,0,1] boundary condition
{for(Iz=1;Iz<=NDZ;Iz++)
{
 for(Iy=1;Iy<=NDY;Iy++)
  {
   for(Ix=1;Ix<=NDX;Ix++)
   {
    E=Iy+NDY*(Ix-1+NDX*(Iz-1));
     N=0;
    for(p=Iz-1;p<=Iz+1;p++)
     {
     for(q=Iy-1;q<=Iy+1;q++)
      {
      for(r=Ix-1;r<=Ix+1;r++)
       {
       N_Iz=p;
       N_Iy=q;
       N_Ix=r;
      
       if(N_Ix==0)
       {N_Ix=NDX;
        break;}
       else if(N_Ix>NDX)
       {N_Ix=1;
        break;}
       if(N_Iy==0)
       {N_Iy=NDY;
       break;}
       else if(N_Iy>NDY)
       {N_Iy=1;
       break;}
       if(N_Iz==0)
       {N_Iz=NDZ;
              }
       else if(N_Iz>NDZ)
       {N_Iz=1;
       }
      else if((N_Ix!=0 || N_Ix<NDX) && (N_Iy!=0 || N_Iy<NDY))
      {
       Neighbor= N_Iy+NDY*(N_Ix-1+NDX*(N_Iz-1));
       Neigh[E][N]=Neighbor;
       N++;
      }
       }
      }
     }
    }
   }
  }
}


if(BC==8)     //[0,0,0] boundary condition
{for(Iz=1;Iz<=NDZ;Iz++)
{
 for(Iy=1;Iy<=NDY;Iy++)
  {
   for(Ix=1;Ix<=NDX;Ix++)
   {
    E=Iy+NDY*(Ix-1+NDX*(Iz-1));
     N=0;
    for(p=Iz-1;p<=Iz+1;p++)
     {
     for(q=Iy-1;q<=Iy+1;q++)
      {
      for(r=Ix-1;r<=Ix+1;r++)
       {
       N_Iz=p;
       N_Iy=q;
       N_Ix=r;
      
       if((N_Ix!=0 || N_Ix<NDX) && (N_Iy!=0 || N_Iy<NDY) && (N_Iz!=0 || N_Iz<NDZ))
      {
       Neighbor= N_Iy+NDY*(N_Ix-1+NDX*(N_Iz-1));
       Neigh[E][N]=Neighbor;
       N++;
      }
      else
      break;
       }
      }
     }
    }
   }
  }
}

//To print the neighbors of 'E'th cell and how many neighbors each cell has.
for(E=1;E<=Totaldomains;E++)
{
 for(N=0;N<27;N++)
 {
   printf("Neighbor of %d is %d \n",E,Neigh[E][N]);
  }
 printf("Number of neighbors is: %d \n",N);
 }

}      
