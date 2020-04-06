#include <coupling1.h>

namespace coupler{

class Part3Data3D{
  public:
    GO  nsurf;    // number of flux surfaces
    GO* versurf; // numbers of vertices on the flux surfaces
//    GO li0,li1,li2,lg1,lg2; // The radial box indexes
    double* xcoords
    GO* xboxinds  //The indexes of all boxes on the radial dimension
    GO* mylk0,myly1,mylk2; // The indexes of box along z dimension
    GO boxsta,boxend,boxcount; // The  indexes of the 2d box

    double** Rcoords;  // The R coordinate of all vertices within the 2d box  
    double** Zcoords;  // The Z coordinate of all vertices within the 2d box 
//    double* density //3D array storing density on each process; 
    double** pzcoords;  // The z coordinates of all points with the 2d box. 

}


void ImportPart3Data3D(Part3Data3D &p3d3d, Part1ParalPar3D  &p1pp3d){
   GO numsurf;
//   GO verosurf;
   if(p1pp3d.mype==0)
     receive_field1D_serial(&numsurf,"../coupling","numsurface",1);
   MPI_Bcast(&numsurf,1,MPI_UNSIGNED_LONG,int root=0, MPI_COMM_WORLD);
   p3d3d.nsurf=numsurf;   
   p3d3d.versurf = new GO[numsurf];
   p3d3d.xcoords = new double[numsurf];
   if(p1pp3d.mype==0){
     receive_field1D_serial(p3d3d.versurf, "../coupling", "vertice_over_surf",numsurf);
     receive_field1D_serial(p3d3d.xcoords,"../coupling", "xcoords_midplane",numsurf);
  }
     MPI_Bcast(p3d3d.versurf,numsurf,MPI_UNSIGNED_LONG,int root=0,MPI_COMM_WORLD);
     MPI_Bcast(p3d3d.x_part3,numsurf,MPI_DOUBLE,int root=0,MPI_COMM_WORLD);
 
   if(preproc==true)
   {
     if(p3d3d.nsurf != p1pp3d.nx0)
     {
       std::cout<<"Error: The number of surface of Part3 doesn't equal to the number vertice of x domain of part1. " \ 
                <<"\n"<<std::endl;
       std::exit;
     }
     GO xinds[]={p1pp3d.li0,p1pp3d.li1,p1pp3d.li2};  
     p3d3d.xboxinds = new GO*[3];
     for(GO i=0;i<3;i++)
     {
       p3d3d.xboxinds[i]=new GO[p1pp3d.npx];
       MPI_Allgather(&xinds[i],1,MPI_UNSIGNED_LONG,p3d3d.xboxinds[i],1,MPI_UNSIGNED_LONG,MPI_COMM_WORLD);
     }       
  } 
}

//when prepro=ture
void DistriPart3zcoords(Part3Data3D &p3d3d, Part1ParalPar3D  &p1pp3d)
{
   GO num=0;
   for(GO i=0;i<p3d3d.nsurf;i++) 
     num+=p3d3d.versurf[i];   
   if(prepro==true)
   {
     double *zcoordall = new double[num]; 
     receive_field1D(zcoordall,"../coupling","all_zcoordinates",num);   
     GO numvert=0, numsurf=0;
     for(GO i=0;i<p1pp3d.mype_x;i++)
     {
       for(GO j=p3d3d.xboxinds[1][i];j<p3d3d.xboxinds[2][i]+1;j++)
       {
         numvert+=p3d3d.versurf[numsurf];
         numsurf+=1; 
       } 
    }
    for(GO i==p3d3d.xboxinds[1][p1pp3d.mype_x];i<p3d3d.xboxinds[2][p1pp3d.mype_x]+1;i++)
    {
       double* zcoords=new double[p3d3d.versurf[numsurf]];
//       numsurf+=1;
       GO numvert1=numvert+p3d3d.versurf[numsurf];
       for(int j=0;j<numvert1;j++)
         zcoords[j]=zcoordall[numvert+j]-pi_;
       GO nstart=minloc(zcoords,p3d3d.versurf[numsurf]);
       reshuffle_nodes(zcoords,nstart,p3d3d.versurf[numsurf]);
           
    }

   }


}

GO  minloc(double* zcoords, const GO n)
  {
    double zmin=std::min(zcoords);
    GO num=0;
    for(GO i=0;i<n;i++) 
     { 
       if(zcoords[i]==zmin) break;
       num=i
     }
     return num;
  }

void reshuffle_nodes(double* zcoords,const GO nstart,const GO vertnum)
{
  double* tmp=new double[vertnum];
  for(GO i=0;i<vertnum-nstart;i++)
    tmp[i]=zcoords[nstart+i];
  for(GO j=vertnum-nstart+1;j<vertnum;j++)
    tmp[j]=zcoords[j-vertnum+nstart-1];
  for(GO k=0;k<vertnum;k++)
    zcoords[k]=tmp[k];
}

void DistributePoints(double* exterarr,int dir,)
{
  

}

void SetDistributeInds(const int dir, li1,li2,double* interarr,double* exterarr)
{
  


}




}
