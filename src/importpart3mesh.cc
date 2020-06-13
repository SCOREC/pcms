#include "importpart3mesh.h"
#include "commpart1.h"
#include "testutilities.h"
#include "sendrecv_impl.h"
#include <cstdlib>
#include <cassert>

namespace coupler{

void Part3Mesh3D::init(const Part1ParalPar3D &p1pp3d,
    const std::string test_dir)
{
   nstart = new LO[p1pp3d.nx0];
   versurf = new LO[p1pp3d.nx0];
   LO numsurf;
   int root=0;
   if(p1pp3d.mype==0){
     if(test_case==TestCase::t0){
       shiftx=-1;
       numsurf=p1pp3d.nx0;
     } else{
        assert(nsurf && cce);
     // Here, numsurf is sent from other parts by Adios routines.  
     }
   }
     if(test_case==TestCase::t0){
        MPI_Bcast(&numsurf,1,MPI_INT,root, MPI_COMM_WORLD);
        nsurf=numsurf; 
      }else{
        MPI_Bcast(&nsurf,1,MPI_INT,root, MPI_COMM_WORLD);
      }
//   if(p1pp3d.mype==0){   // please keep this commented loop.
   if(test_case==TestCase::t0){
      versurfpart3 = new LO[nsurf];
      xcoords = new double[nsurf];

      assert(!test_dir.empty());
      std::string fname=test_dir+"versurf.nml";
      InputfromFile(versurfpart3,numsurf,fname);
      fname=test_dir+"xcoords.nml";
      InputfromFile(xcoords,nsurf,fname);
// please keep the following two commented lines.After determining the communnicator, they will be removed.
//       MPI_Bcast(versurfpart3,nsurf,MPI_INT,root,MPI_COMM_WORLD);
//       MPI_Bcast(xcoords,nsurf,MPI_DOUBLE,root,MPI_COMM_WORLD);

      cce_first_surface=0;
      cce_last_surface=nsurf-1;
      cce_first_node=1;
      cce_last_node=0;
      for(LO i=0;i<nsurf;i++)
	cce_last_node+=(GO)versurfpart3[i];
      cce_node_number = cce_last_node-cce_first_node+1;

    }else {

      cce_first_surface=(LO)cce[0];
      cce_last_surface=(LO)cce[1];
      cce_first_node=cce[2];
      cce_last_node=cce[3];
      cce_node_number = cce_last_node-cce_first_node+1;
      shiftx = cce_first_surface-1; 
      assert(versurfpart3);
      assert(xcoords);
    }
//   }


   if(preproc==true){
     JugeFirstSurfaceMatch(p1pp3d.xcoords[0]); // Judge whether the first surface of part1 
                                               // equals cce_first_surface of part3 
     for(LO i=0;i<p1pp3d.nx0; i++)
       versurf[i]=versurfpart3[i+shiftx+1];     

     totnode=0;
     for(LO i=0;i<nsurf;i++)
       totnode+=(GO)versurfpart3[i];
  
     activenode=0;
     std::cerr<<"nx0: "<<p1pp3d.nx0<<"\n";
     for(LO i=0;i<p1pp3d.nx0;i++){
      MPI_Barrier(MPI_COMM_WORLD);
     //std::cerr<<"i: "<<i<<" versurf[i]: "<<versurf[i]<<"\n";
       activenode+=(GO)versurf[i];}
    std::cerr<<"activenode: "<<activenode<<" cce_node_number "<<cce_node_number<<"\n";
     if(activenode!=cce_node_number){
       std::cout<<"ERROR: The activenode number of part1 doesn't equal to cce_node_number for part3."<<'\n';
       std::exit(EXIT_FAILURE);
     }
      std::cerr<<"0.0"<<"\n"; 
     li0=p1pp3d.li0;
     li1=p1pp3d.li1;
     li2=p1pp3d.li2;

      std::cerr<<"0.1"<<"\n"; 
     BlockIndexes(p1pp3d.comm_x,p1pp3d.mype_x,p1pp3d.npx); 
 
      std::cerr<<"0.2"<<"\n"; 
     LO xinds[3]={p1pp3d.li0,p1pp3d.li1,p1pp3d.li2}; 
     xboxinds = new LO*[p1pp3d.npx]; 
     for(LO i=0;i<p1pp3d.npx;i++){
       xboxinds[i]=new LO[3];
     }
      std::cerr<<"0.3"<<"\n"; 

     LO* buffer=new LO[3*p1pp3d.npx];
     MPI_Allgather(xinds,3,MPI_INT,buffer,3,MPI_INT,p1pp3d.comm_x);
     for(LO i=0;i<p1pp3d.npx;i++){
       for(LO j=0;j<3;j++)
         xboxinds[i][j]=buffer[i*3+j];
     }
      std::cerr<<"0.4"<<"\n"; 
     delete[] buffer;
     if(test_case==TestCase::t0){
       if(p1pp3d.mype_x==0){
	 for(LO k=0;k<3;k++){
      	     std::cout<<"xboxinds[1]["<<k<<"]="<<xboxinds[1][k]<<'\n';
         }
       }
     }
      std::cerr<<"0.5"<<"\n"; 
     lj0=p1pp3d.lj0*2; 
     // FIXME mylk0 is undersized for circular case;
     // it is written in DistributePoints and
     // read in DistriPart3zcoords
     mylk0=new LO[li0]; 
     mylk1=new LO[li0];
     mylk2=new LO[li0];      
      std::cerr<< p1pp3d.mype << " 0.6"<<"\n"; 
     DistriPart3zcoords(p1pp3d, test_dir);
      std::cerr<< p1pp3d.mype << " 0.7"<<"\n"; 
  } 
}

void Part3Mesh3D::BlockIndexes(const MPI_Comm comm_x,const LO mype_x,const LO npx)
{
  GO* inds = new GO[npx]; 
  blockcount=0;
  for(LO i=0;i<li0;i++)
    blockcount+=(GO)versurf[li1+i];
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
  MPI_Barrier(MPI_COMM_WORLD);
  if(!p1pp3d.mype)  std::cerr<<"0.60"<<"\n"; 
  if(preproc==true){
    if(test_case==TestCase::t0){
      zcoordall = new double[activenode];
      InitzcoordsInCoupler(zcoordall,versurf,p1pp3d.nx0);
    }else{
      MPI_Barrier(MPI_COMM_WORLD);
      if(!p1pp3d.mype)  std::cerr<<"0.61"<<"\n"; 
      assert(zcoordall); // the number of elements is activenode
    }
    LO numvert=0, numsurf=0;
    MPI_Barrier(MPI_COMM_WORLD);
    if(!p1pp3d.mype)  std::cerr<<"0.615"<<"\n"; 
    for(LO i=0;i<p1pp3d.mype_x;i++){
      if(!p1pp3d.mype)  std::cerr<<"0.619"<<"\n"; 
      for(LO j=xboxinds[i][1];j<xboxinds[i][2]+1;j++){
        std::cerr<<p1pp3d.mype<<" 0.62,numsurf: "<<numsurf<<" versurf[numsurf] "<<versurf[numsurf]<<"\n";
	numvert+=versurf[numsurf];
	numsurf+=1; 
      } 
    }
    MPI_Barrier(MPI_COMM_WORLD);
    std::cerr<<"0.63"<<"\n"; 
    LO index1=xboxinds[p1pp3d.mype_x][1];
    LO index2=xboxinds[p1pp3d.mype_x][2];
    LO index0=xboxinds[p1pp3d.mype_x][0];
    pzcoords = new double*[index0]; 
    MPI_Barrier(MPI_COMM_WORLD);
    std::cerr<<p1pp3d.mype << " 0.64, numsurf: "<<numsurf<<" numvert: "<<numvert<<"\n"; 
    double* zcoords;
    std::cerr<< p1pp3d.mype << " loop range " << index1 << " : " << index2+1 << "\n";
    MPI_Barrier(MPI_COMM_WORLD);
    for(LO i= index1;i<index2+1;i++)
    {
      std::cerr<< p1pp3d.mype << " 0.641"<<" versurf[numsurf] "<<versurf[numsurf]<<"\n"; 
      zcoords=new double[versurf[numsurf]];  
      std::cerr<< p1pp3d.mype << " 0.642"<<"\n"; 
      for(LO j=0;j<versurf[numsurf];j++){
      std::cerr<< p1pp3d.mype << " 0.643, j "<< j <<" numvert+j "<< numvert+j <<" zcoordall[numvert+j] "<<zcoordall[numvert+j]<<"\n"; 
      std::cerr<< p1pp3d.mype << " 0.643, numvert+j "<< numvert+j <<"\n"; 
	 zcoords[j]=zcoordall[numvert+j]-cplPI;
      }
      std::cerr<< p1pp3d.mype << " 0.65"<<"\n"; 
      if(test_case==TestCase::t0){
        assert(!test_dir.empty());
        std::string fname=test_dir+std::to_string(i)+"_zcoords.txt";
        OutputtoFile(zcoords,versurf[numsurf],fname);
      }
      std::cerr<<p1pp3d.mype << " idx " << i << " 0.66"<<"\n"; 
      nstart[i] = minloc(zcoords,versurf[numsurf]);
/*
if(p1pp3d.mype_x==1){
  std::cout<<"versurf["<<numsurf<<"]="<<versurf[numsurf]<<'\n';
  std::cout<<"li="<<i<<" "<<"nstart="<<nstart<<'\n';
}
*/
       reshuffleforward(zcoords,nstart[i],versurf[numsurf]);
       std::cerr<<p1pp3d.mype << " 0.67,numsurf: "<<numsurf<<" versurf[numsurf] "<<versurf[numsurf]<<"\n"; 
       DistributePoints(zcoords,index1,i,p1pp3d.pzcoords,p1pp3d);
       std::cerr<<p1pp3d.mype << " 0.68 p1.nz0 " << p1pp3d.nz0 << "mylk0[i-index1] " << mylk0[i-index1] << "\n";
       pzcoords[i-index1]= new double[mylk0[i-index1]]; // FIXME li is ~0:90 for circular case, this read is out of bounds
       std::cerr<<p1pp3d.mype << " 0.69"<<"\n"; 
       for(LO k=0;k<mylk0[i-index1];k++){ // FIXME li is ~0:90 for circular case, this read is out of bounds
	 pzcoords[i-index1][k]= zcoords[mylk1[i-index1]+k];
       }
       std::cerr<<p1pp3d.mype << " 0.691"<<"\n"; 
       numvert+=versurf[numsurf];
       numsurf+=1;       
       std::cerr<<p1pp3d.mype << " 0.692"<<"\n"; 
       delete[] zcoords; 
       std::cerr<<p1pp3d.mype << " idx " << i << " 0.693"<<"\n"; 
    }
    std::cerr<<p1pp3d.mype << " 0.6931"<<"\n"; 
  }  
  std::cerr<<p1pp3d.mype << " 0.694"<<"\n"; 
}

void Part3Mesh3D::JugeFirstSurfaceMatch(double xp1)
{
  double* tmp = new double[nsurf];
  for(LO i=0;i<nsurf; i++){
    tmp[i]=abs(xcoords[i]-(xp1));
  }
  if(test_case==TestCase::t0){
    shiftx = minloc(tmp,nsurf)-1;
  }else{};
  std::cout<<"cce_first_surface "<<cce_first_surface<<'\n';
  if(shiftx+1!=cce_first_surface){
    std::cout<<"shiftx="<<shiftx<<'\n';
    std::cout<<"ERROR: The first surface of part1 doesn't match cce_first_surface."<<'\n';
    exit(1);
  }
  delete[] tmp;
}

LO  minloc(const double* array, const LO n)
{
    double zmin=minimalvalue(array, n);
    LO num=0;
    for(LO i=0;i<n;i++){ 
       if(array[i]==zmin) break;
       num=i;
     }
     return num;
 }

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
    std::cerr << p1pp3d.mype << " li-gstart " << li-gstart << "\n";
    mylk0[li-gstart]=i2-i1+1; // li is ~0:90 for circular case, this write is out of bounds
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
