#include "importpart3mesh.h"
#include "commpart1.h"
#include "testutilities.h"
#include "sendrecv_impl.h"
#include <cassert>

namespace coupler{

void Part3Mesh3D::init(const Part1ParalPar3D &p1pp3d,
    const std::string test_dir)
{
   LO numsurf;
   int root=0;
   if(p1pp3d.mype==0){
     if(test_case==TestCase::t0){
       numsurf=p1pp3d.nx0;
     } else{
//     receive_field1D_serial(&numsurf,"../coupling","numsurface",1); 
     }
   }
   MPI_Bcast(&numsurf,1,MPI_INT,root, MPI_COMM_WORLD);
   nsurf=numsurf; 
   nstart = new LO[numsurf];  
   versurf = new LO[numsurf];
   
   xcoords = new double[numsurf];
   if(p1pp3d.mype==0){
     if(test_case==TestCase::t0){
       assert(!test_dir.empty());
       std::string fname=test_dir+"versurf.nml";
       InputfromFile(versurf,numsurf,fname);
       fname=test_dir+"xcoords.nml";
       InputfromFile(xcoords,numsurf,fname);
     }else {
//     receive_field1D_serial(versurf, "../coupling", "vertice_over_surf",numsurf);
//     receive_field1D_serial(xcoords,"../c
//     upling", "xcoords_midplane",numsurf);
     }
   }
   MPI_Bcast(versurf,numsurf,MPI_INT,root,MPI_COMM_WORLD);
   MPI_Bcast(xcoords,numsurf,MPI_DOUBLE,root,MPI_COMM_WORLD);

   totnode=0;
   for(LO i=0;i<nsurf;i++)
     totnode+=(GO)versurf[i];

   if(preproc==true){
     if(nsurf != p1pp3d.nx0)
     {std::cout<<"Error: The number of surface of Part3 doesn't equal to the number vertice of x domain of part1. "
               <<"\n"<<std::endl;
       std::exit(EXIT_FAILURE);
     }
     li0=p1pp3d.li0;
     li1=p1pp3d.li1;
     li2=p1pp3d.li2;

     BlockIndexes(p1pp3d.comm_x,p1pp3d.mype_x,p1pp3d.npx); 
 
     LO xinds[3]={p1pp3d.li0,p1pp3d.li1,p1pp3d.li2}; 
     xboxinds = new LO*[p1pp3d.npx]; 
     for(LO i=0;i<p1pp3d.npx;i++){
       xboxinds[i]=new LO[3];
     }
     LO* buffer=new LO[3*p1pp3d.npx];
     MPI_Allgather(xinds,3,MPI_INT,buffer,3,MPI_INT,p1pp3d.comm_x);
     for(LO i=0;i<p1pp3d.npx;i++){
       for(LO j=0;j<3;j++)
         xboxinds[i][j]=buffer[i*3+j];
     }
     delete[] buffer;
     if(test_case==TestCase::t0){
       if(p1pp3d.mype_x==0){
	 for(LO k=0;k<3;k++){
      	     std::cout<<"xboxinds[1]["<<k<<"]="<<xboxinds[1][k]<<'\n';
         }
       }
     }
     lj0=p1pp3d.lj0*2; 
     mylk0=new LO[3];
     mylk1=new LO[li0];
     mylk2=new LO[li0];      
     DistriPart3zcoords(p1pp3d, test_dir);
  } 
}

void Part3Mesh3D::BlockIndexes(const MPI_Comm comm_x,const LO mype_x,const LO npx)
{
  GO* inds = new GO[npx]; 
  blockcount=0;
  for(LO i=0;i<li0;i++)
    blockcount+=blockcount+(GO)versurf[li1+i];
  MPI_Datatype mpitype;
  mpitype = getMpiType(GO());
  MPI_Allgather(MPI_IN_PLACE,1,mpitype,inds,1,mpitype,comm_x);
  blockstart=0;
  for(LO i=0;i<mype_x;i++) 
     blockstart+=inds[i];
  blockend=blockstart+blockcount-1;
  delete[] inds;
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

//when prepro=ture
void Part3Mesh3D::DistriPart3zcoords(const Part1ParalPar3D &p1pp3d,
    const std::string test_dir)
{
  if(preproc==true){
    GO num=0;
    for(LO i=0;i<nsurf;i++) 
       num+=(GO)(versurf[i]); 
    double *zcoordall = new double[num];
    if(test_case==TestCase::t0){
      InitzcoordsInCoupler(zcoordall,versurf,nsurf);
    }else{
//     receive_field1D(zcoordall,"../coupling","all_zcoordinates",num);   
    }
    LO numvert=0, numsurf=0;
    for(LO i=0;i<p1pp3d.mype_x;i++){
      for(LO j=xboxinds[i][1];j<xboxinds[i][2]+1;j++){
	numvert+=versurf[numsurf];
	numsurf+=1; 
      } 
    }

    LO index1=xboxinds[p1pp3d.mype_x][1];
    LO index2=xboxinds[p1pp3d.mype_x][2];
    LO index0=xboxinds[p1pp3d.mype_x][0];
    pzcoords = new double*[index0]; 
    double* zcoords;
    for(LO i= index1;i<index2+1;i++)
    {
      zcoords=new double[versurf[numsurf]];  
      for(LO j=0;j<versurf[numsurf];j++){
	 zcoords[j]=zcoordall[numvert+j]-cplPI;
      }
      if(test_case==TestCase::t0){
        assert(!test_dir.empty());
        std::string fname=test_dir+std::to_string(i)+"_zcoords.txt";
       OutputtoFile(zcoords,versurf[numsurf],fname);
      }
      nstart[i] = minloc(zcoords,versurf[numsurf]);
/*
if(p1pp3d.mype_x==1){
  std::cout<<"versurf["<<numsurf<<"]="<<versurf[numsurf]<<'\n';
  std::cout<<"li="<<i<<" "<<"nstart="<<nstart<<'\n';
}
*/
       reshuffleforward(zcoords,nstart[i],versurf[numsurf]);
       DistributePoints(zcoords,index1,i,p1pp3d.pzcoords,p1pp3d);
       pzcoords[i-index1]= new double[mylk0[i-index1]];
       for(LO k=0;k<mylk0[i-index1];k++){
	 pzcoords[i-index1][k]= zcoords[mylk1[i-index1]+k];
       }
       numvert+=versurf[numsurf];
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
/*
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
*/

//// notice: be carefull with extra_zero case.
// CWS - one of the classes must be read only.... 
void Part3Mesh3D::DistributePoints(const double* exterarr, const LO gstart,LO li, 
                  const double* interarr, const Part1ParalPar3D  &p1pp3d)
{
  if(preproc==true){
    LO nstart;
    double* tmp=new double[versurf[li]];
    for(LO i=0;i<versurf[li];i++)
      tmp[i]=abs(exterarr[i]-interarr[p1pp3d.lk1]);
    nstart=minloc(tmp,versurf[li]);
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
      if(i2>=versurf[li]-1){
        break;
      }
      if(exterarr[i2+1]<internal_ub){
        i2+=1;
      }
      else{
        inside=false;
      }
    }
    mylk1[li-gstart]=i1;
    mylk2[li-gstart]=i2;
    mylk0[li-gstart]=i2-i1+1;
    if(test_case==TestCase::t0){   
      std::cout<<"rank="<<p1pp3d.mype<<" "<<li-gstart<<'\n';
      std::cout<<"mylk k="<<mylk0[li-gstart]<<" "<<mylk1[li-gstart]
      <<" "<<mylk2[li-gstart]<<" "<<'\n'; 
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


}
