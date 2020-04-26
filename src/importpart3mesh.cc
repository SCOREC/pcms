#include "importpart3mesh.h"

namespace coupler{

void ImportPart3Mesh3D(Part3Mesh3D &p3m3d, Part1ParalPar3D  &p1pp3d)
{
   LO numsurf;
   int root=0;
   if(p1pp3d.mype==0){
     if(test_case==0){
       numsurf=p1pp3d.nx0;
     } else{
//     receive_field1D_serial(&numsurf,"../coupling","numsurface",1); 
     }
   }
   MPI_Bcast(&numsurf,1,MPI_INT,root, MPI_COMM_WORLD);
   p3m3d.nsurf=numsurf;   
   p3m3d.versurf = new LO[numsurf];
   p3m3d.xcoords = new double[numsurf];
   if(p1pp3d.mype==0){
     if(test_case==0){
       std::string fname=test_dir+"versuf.nml";
       InputfromFile(p3m3d.versurf,numsurf,fname);
       fname=test_dir+"xcoords.nml";
       InputfromFile(p3m3d.xcoords,numsurf,fname);
     }else {
//     receive_field1D_serial(p3m3d.versurf, "../coupling", "vertice_over_surf",numsurf);
//     receive_field1D_serial(p3m3d.xcoords,"../c
//     upling", "xcoords_midplane",numsurf);
     }
   }
     MPI_Bcast(p3m3d.versurf,numsurf,MPI_UNSIGNED_LONG,root,MPI_COMM_WORLD);
     MPI_Bcast(p3m3d.xcoords,numsurf,MPI_DOUBLE,root,MPI_COMM_WORLD);
 
   if(preproc==true){
     if(p3m3d.nsurf != p1pp3d.nx0)
     {std::cout<<"Error: The number of surface of Part3 doesn't equal to the number vertice of x domain of part1. " \ 
               <<"\n"<<std::endl;
       std::exit(EXIT_FAILURE);
     }
     p3m3d.li0=p1pp3d.li0;
     p3m3d.li1=p1pp3d.li1;
     p3m3d.li2=p1pp3d.li2;
     LO xinds[]={p1pp3d.li0,p1pp3d.li1,p1pp3d.li2};  
     p3m3d.xboxinds = new LO*[3];
     for(LO i=0;i<3;i++){
       p3m3d.xboxinds[i]=new LO[p1pp3d.npx];
       MPI_Allgather(&xinds[i],1,MPI_UNSIGNED_LONG,p3m3d.xboxinds[i],1,MPI_UNSIGNED_LONG,MPI_COMM_WORLD);
     }
     p3m3d.lj0=p1pp3d.lj0*2; 
     p3m3d.mylk0=new LO[p3m3d.li0];
     p3m3d.mylk1=new LO[p3m3d.li0];
     p3m3d.mylk2=new LO[p3m3d.li0];      
     DistriPart3zcoords(p3m3d,p1pp3d);
  } 
}

//when prepro=ture
void DistriPart3zcoords(Part3Mesh3D &p3m3d, Part1ParalPar3D  &p1pp3d)
{
  if(preproc==true){
    LO num=0;
    for(LO i=0;i<p3m3d.nsurf;i++) 
       num+=p3m3d.versurf[i];   
    if(preproc==true){
      double *zcoordall = new double[num];
      if(test_case==0){
  //      InputfromFile(zcoordall,num,"zcoordall.rtf");
	InitzcoordsInCoupler(zcoordall,num);
      }else{ 
  //     receive_field1D(zcoordall,"../coupling","all_zcoordinates",num);   
      }
      LO numvert=0, numsurf=0;
      for(LO i=0;i<p1pp3d.mype_x;i++){
	for(LO j=p3m3d.xboxinds[1][i];j<p3m3d.xboxinds[2][i]+1;j++){
	  numvert+=p3m3d.versurf[numsurf];
	  numsurf+=1; 
	} 
      }

      LO index1=p3m3d.xboxinds[1][p1pp3d.mype_x];
      LO index2=p3m3d.xboxinds[2][p1pp3d.mype_x];
      LO index0=p3m3d.xboxinds[0][p1pp3d.mype_x];
      p3m3d.pzcoords = new double*[index0]; 
      double* zcoords;
      for(LO i= index1;i<index2+1;i++)
      {
	zcoords=new double[p3m3d.versurf[numsurf]];  
	for(int j=0;j<p3m3d.versurf[numsurf];j++){
	   zcoords[j]=zcoordall[numvert+j]-cplPI;
	}
  //       numvert+=p3m3d.versurf[numsurf];
	 LO nstart=minloc(zcoords,p3m3d.versurf[numsurf]);
	 reshuffle_nodes(zcoords,nstart,p3m3d.versurf[numsurf]);
	 DistributePoints(zcoords,index1,i,p1pp3d.pzcoords,p3m3d,p1pp3d);
	 p3m3d.pzcoords[i-index1]= new double[p3m3d.mylk0[i-index1]];
	 for(LO k=0;k<p3m3d.mylk0[i-index1];k++){
	   p3m3d.pzcoords[i-index1][k]= zcoords[p3m3d.mylk1[i-index1]+k];
	 }
	 numvert+=p3m3d.versurf[numsurf];
	 numsurf+=1;       
	 delete[] zcoords; 
      }
   }
 }  
}

LO  minloc(const double* zcoords, const LO n)
{
    double zmin=minimalvalue(zcoords, n);
    LO num=0;
    for(LO i=0;i<n;i++){ 
       if(zcoords[i]==zmin) break;
       num=i;
     }
     return num;
 }

void reshuffle_nodes(double* zcoords,const LO nstart,const LO vertnum)
{
  double* tmp=new double[vertnum];
  for(LO i=0;i<vertnum-nstart;i++)
    tmp[i]=zcoords[nstart+i];
  for(LO j=vertnum-nstart+1;j<vertnum;j++)
    tmp[j]=zcoords[j-vertnum+nstart-1];
  for(LO k=0;k<vertnum;k++)
    zcoords[k]=tmp[k];
}

//// notice: be carefull with extra_zero case.
void DistributePoints(double* exterarr,LO gstart,LO li, double* interarr,Part3Mesh3D &p3m3d, Part1ParalPar3D  &p1pp3d)
{
  if(preproc==true){
    LO nstart;
    double* tmp=new double[p3m3d.versurf[li]];
    for(LO i=0;i<p3m3d.versurf[li];i++)
      tmp[i]=abs(exterarr[i]-interarr[p1pp3d.lk1]);
    nstart=minloc(tmp,p3m3d.versurf[li]);
    //nstart must be in my domain or will duplicate
    if(exterarr[nstart]<interarr[p1pp3d.lk1])
      nstart+=1;
    LO i1=nstart;
    LO i2=nstart;
    double internal_ub;
    if(p1pp3d.lk2==p1pp3d.nz0-1){
      internal_ub=cplPI;
    }
    else{
      internal_ub=p1pp3d.pzcoords[p1pp3d.lk2+1];
    }
    bool inside = true;
    while(inside){
      if(i2+1>p3m3d.versurf[li]){
        break;
      }
      if(exterarr[i2+1]<internal_ub){
        i2+=1;
      }
      else{
        inside=false;
      }
    }
    p3m3d.mylk1[li-gstart]=i1;
    p3m3d.mylk2[li-gstart]=i2;
    p3m3d.mylk0[li-gstart]=i2-i1+1;
    }
}

 double minimalvalue(const double* array, const LO n)
{
    double tmp=array[0];
    for(LO i=1;i<n;i++){
      if(tmp>array[i]){
        tmp=array[i];
      }
    }
    return tmp;      
}
 
void InitzcoordsInCoupler(double* zcoords,LO num)
{
  double shift=0.1;
  double delta=2.0*cplPI/(double)num;
  for(LO i=0;i<num-1;i++){
    zcoords[i]=double(i)*delta+shift;
  }
} 


}


