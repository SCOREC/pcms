#include <coupling1.h>

namespace coupler{

class Part3Mesh3D{
  public:
    GO  nsurf;    // number of flux surfaces
    GO* versurf; // numbers of vertice on the flux surfaces
//    GO li0,li1,li2,lg1,lg2; // The radial box indexes
    double* xcoords;
    GO  li0,li1,li2;
    GO** xboxinds;  //The indexes of all boxes on the radial dimension
    GO lj0;
    GO* mylk0,mylk1,mylk2; // The indexes of box along z dimension
    GO boxstar,boxend,boxcount; // The  indexes of the 2d box

    double** Rcoords;  // The R coordinate of all vertices within the 2d box  
    double** Zcoords;  // The Z coordinate of all vertices within the 2d box 
    double** pzcoords;  // The z coordinates of all points with the 2d box. 

}


void ImportPart3Mesh3D(Part3Mesh3D &p3m3d, Part1ParalPar3D  &p1pp3d){
   GO numsurf;
//   GO verosurf;
   if(p1pp3d.mype==0)
     receive_field1D_serial(&numsurf,"../coupling","numsurface",1);
   MPI_Bcast(&numsurf,1,MPI_UNSIGNED_LONG,int root=0, MPI_COMM_WORLD);
   p3m3d.nsurf=numsurf;   
   p3m3d.versurf = new GO[numsurf];
   p3m3d.xcoords = new double[numsurf];
   if(p1pp3d.mype==0){
     receive_field1D_serial(p3m3d.versurf, "../coupling", "vertice_over_surf",numsurf);
     receive_field1D_serial(p3m3d.xcoords,"../coupling", "xcoords_midplane",numsurf);
  }
     MPI_Bcast(p3m3d.versurf,numsurf,MPI_UNSIGNED_LONG,int root=0,MPI_COMM_WORLD);
     MPI_Bcast(p3m3d.x_part3,numsurf,MPI_DOUBLE,int root=0,MPI_COMM_WORLD);
 
   if(preproc==true){
     if(p3m3d.nsurf != p1pp3d.nx0)
     {std::cout<<"Error: The number of surface of Part3 doesn't equal to the number vertice of x domain of part1. " \ 
               <<"\n"<<std::endl;
       std::exit;
     }
     p3m3d.li0=p1pp3d.li0;
     p3m3d.li1=p1pp3d.li1;
     p3m3d.li2=p1pp3d.li2;
     GO xinds[]={p1pp3d.li0,p1pp3d.li1,p1pp3d.li2};  
     p3m3d.xboxinds = new GO*[3];
     for(GO i=0;i<3;i++){
       p3m3d.xboxinds[i]=new GO[p1pp3d.npx];
       MPI_Allgather(&xinds[i],1,MPI_UNSIGNED_LONG,p3m3d.xboxinds[i],1,MPI_UNSIGNED_LONG,MPI_COMM_WORLD);
     }       
   DistriPart3zcoords(&p3m3d,&p1pp3d);
  } 
}

//when prepro=ture
void DistriPart3zcoords(Part3Mesh3D &p3m3d, Part1ParalPar3D  &p1pp3d){
   GO num=0;
   for(GO i=0;i<p3m3d.nsurf;i++) 
     num+=p3m3d.versurf[i];   
   if(prepro==true){
     double *zcoordall = new double[num]; 
     receive_field1D(zcoordall,"../coupling","all_zcoordinates",num);   
     GO numvert=0, numsurf=0;
     for(GO i=0;i<p1pp3d.mype_x;i++){
       for(GO j=p3m3d.xboxinds[1][i];j<p3m3d.xboxinds[2][i]+1;j++){
         numvert+=p3m3d.versurf[numsurf];
         numsurf+=1; 
       } 
    }

    p3m3d.pzcoords = new double*[p3m3d.xboxinds[0][p2pp3d.mype_x]];
    GO index=p3m3d.xboxinds[1][p1pp3d.mype_x];
    for(GO i== index;i<p3m3d.xboxinds[2][p1pp3d.mype_x]+1;i++)
    {
       double* zcoords=new double[p3m3d.versurf[numsurf]];
//       numsurf+=1;
       GO numvert1=numvert+p3m3d.versurf[numsurf];
       for(int j=0;j<numvert1;j++)
         zcoords[j]=zcoordall[numvert+j]-pi_;
       GO nstart=minloc(zcoords,p3m3d.versurf[numsurf]);
       reshuffle_nodes(zcoords,nstart,p3m3d.versurf[numsurf]);
       DistributePoints(zcoords,index,i,p1pp3d.pzcoords,&p3m3d,&p1pp3d)
       p3m3d.pzcoords[i-index]= new double*[p3m3d.mylk0[i-index]];
       for(GO k=0;k<p3m3d.mylk0[i-index];k++)
         p3m3d.pzcoords[i-index][k]=zcoords[p3m3d.mylk1[i-index]+k];
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

//// notice: be carefull with extra_zero case.
void DistributePoints(double* exterarr,GO gstart,GO li, double* interarr,Part3Mesh3D &p3m3d, Part1ParalPar3D  &p1pp3d){
  if(prepro==true){
    GO nstart;
    double* tmp=new double[p3m3d.versurf[li]]
    for(GO i=0;i<p3m3d.versurf[li];i++)
      tmp[i]=abs(exterarr[i]-interarr[p1pp3d.lk1]);
    nstart=minloc(tmp,p1pp3d.li0);
    //nstart must be in my domain or will duplicate
    if(exterarr[nstart]<interarr[p1pp3d.lk1])
      nstart+=1;
    GO i1=nstart;
    GO i2=nstart;
    double interal_ub;
    if(p1pp3d.lk2==p1pp3d.nz0-1){
      internal_ub=pi_;
    }
    else{
      internal_ub=p1pp3d.pzcoords[p1pp3d.lk2+1];
    }
    bool inside = true;
    while(bool){
      if(i2+1>p1pp3d.versurf[li]){
        break;
      }
      if(exterarr[i2+1]<interal_ub){
        i2+=1;
      }
      else{
        inside=false;
      }
    }
    p3m3d.mylk1[li-gstart]=i1;
    p3m3d.mylk2[li-gstart]=i2;
    p3m3d.mylk0[li-gstart]=i2-i2+1;
    }
}

}
