#include "importpart3mesh.h"
#include "IOutilities.h"

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
       std::string fname=test_dir+"versurf.nml";
       InputfromFile(p3m3d.versurf,numsurf,fname);
       fname=test_dir+"xcoords.nml";
       InputfromFile(p3m3d.xcoords,numsurf,fname);
     }else {
//     receive_field1D_serial(p3m3d.versurf, "../coupling", "vertice_over_surf",numsurf);
//     receive_field1D_serial(p3m3d.xcoords,"../c
//     upling", "xcoords_midplane",numsurf);
     }
   }
     MPI_Bcast(p3m3d.versurf,numsurf,MPI_INT,root,MPI_COMM_WORLD);
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
     LO xinds[3]={p1pp3d.li0,p1pp3d.li1,p1pp3d.li2};  
     p3m3d.xboxinds = new LO*[p1pp3d.npx]; 
     for(LO i=0;i<p1pp3d.npx;i++){
       p3m3d.xboxinds[i]=new LO[3];
     }
     LO* buffer=new LO[3*p1pp3d.npx];
     MPI_Allgather(xinds,3,MPI_INT,buffer,3,MPI_INT,p1pp3d.comm_x);
     for(LO i=0;i<p1pp3d.npx;i++){
       for(LO j=0;j<3;j++)
         p3m3d.xboxinds[i][j]=buffer[i*3+j];
     }
     delete[] buffer;
     if(test_case==0){
       if(p1pp3d.mype_x==0){
	 for(LO k=0;k<3;k++)
	 std::cout<<"k"<<" "<<p3m3d.xboxinds[1][k]<<'\n';
       }
     }
     p3m3d.lj0=p1pp3d.lj0*2; 
     p3m3d.mylk0=new LO[3];
     p3m3d.mylk1=new LO[p3m3d.li0];
     p3m3d.mylk2=new LO[p3m3d.li0];      
     DistriPart3zcoords(p3m3d,p1pp3d);
  } 
}

//when prepro=ture
void DistriPart3zcoords(Part3Mesh3D &p3m3d, Part1ParalPar3D  &p1pp3d)
{
  if(preproc==true){
    GO num=0;
    for(LO i=0;i<p3m3d.nsurf;i++) 
       num+=(GO)(p3m3d.versurf[i]); 
    double *zcoordall = new double[num];
    if(test_case==0){
      InitzcoordsInCoupler(zcoordall,p3m3d.versurf,p3m3d.nsurf);
    }else{
//     receive_field1D(zcoordall,"../coupling","all_zcoordinates",num);   
    }
    LO numvert=0, numsurf=0;
    for(LO i=0;i<p1pp3d.mype_x;i++){
      for(LO j=p3m3d.xboxinds[i][1];j<p3m3d.xboxinds[i][2]+1;j++){
	numvert+=p3m3d.versurf[numsurf];
	numsurf+=1; 
      } 
    }

    LO index1=p3m3d.xboxinds[p1pp3d.mype_x][1];
    LO index2=p3m3d.xboxinds[p1pp3d.mype_x][2];
    LO index0=p3m3d.xboxinds[p1pp3d.mype_x][0];
    p3m3d.pzcoords = new double*[index0]; 
    double* zcoords;
    for(LO i= index1;i<index2+1;i++)
    {
      zcoords=new double[p3m3d.versurf[numsurf]];  
      for(LO j=0;j<p3m3d.versurf[numsurf];j++){
	 zcoords[j]=zcoordall[numvert+j]-cplPI;
      }
      if(test_case==0){
        std::string fname=test_dir+std::to_string(i)+"_zcoords.txt";
       OutputtoFile(zcoords,p3m3d.versurf[numsurf],fname);
      }
       LO nstart=minloc(zcoords,p3m3d.versurf[numsurf]);
/*
if(p1pp3d.mype_x==1){
  std::cout<<"p3m3d.versurf["<<numsurf<<"]="<<p3m3d.versurf[numsurf]<<'\n';
  std::cout<<"li="<<i<<" "<<"nstart="<<nstart<<'\n';
}
*/
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
    if(test_case==0){   
      std::cout<<"rank="<<p1pp3d.mype<<" "<<li-gstart<<'\n';
      std::cout<<"mylk k="<<p3m3d.mylk0[li-gstart]<<" "<<p3m3d.mylk1[li-gstart] \
      <<" "<<p3m3d.mylk2[li-gstart]<<" "<<'\n'; 
    }
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
 
void InitzcoordsInCoupler(double* zcoords,LO* versurf,LO nsurf)
{
  double shift=0.1;
  GO num=0;
  for(LO i=0;i<nsurf;i++){
    double delta=2.0*cplPI/(double)(versurf[i]);
    for(LO j=0;j<versurf[i];j++){
      zcoords[num]=(double)(j)*delta+shift;
      num++;
    }
  } 
}



}
