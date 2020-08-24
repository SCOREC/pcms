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
   PERFSTUBS_START_STRING(__func__);
   nstart = new LO[p1pp3d.nx0];
   LO numsurf;
   int root=0;
   if(p1pp3d.mype==0){
     if(test_case==TestCase::t0){
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
   if(test_case==TestCase::t0){
      versurf = new LO[nsurf];
      xcoords = new double[nsurf];

      assert(!test_dir.empty());
      std::string fname=test_dir+"versurf.nml";
      InputfromFile(versurf,numsurf,fname);
      fname=test_dir+"xcoords.nml";
      InputfromFile(xcoords,nsurf,fname);

      cce_first_surface=0;
      cce_last_surface=nsurf-1;
      cce_first_node=1;
      cce_last_node=0;
      shiftx=-1; 
      for(LO i=0;i<nsurf;i++)
	cce_last_node+=(GO)versurf[i];
      cce_node_number = cce_last_node-cce_first_node+1;

    }else {

      cce_first_surface=(LO)cce[0];
      cce_last_surface=(LO)cce[1];
      cce_first_node=cce[2];
      cce_last_node=cce[3];
      cce_node_number = cce_last_node-cce_first_node+1;
      shiftx = cce_first_surface-1; 
      assert(versurf);
      assert(xcoords);

    }


   if(preproc==true){
     JugeFirstSurfaceMatch(p1pp3d.xcoords[0]); // Judge whether the first surface of part1 

     activenodes=0;
     for(LO i=0;i<p1pp3d.nx0;i++){
       activenodes+=(GO)versurf[i];
     }

     if(p1pp3d.mype==0){
       std::cout<<"activenodes, cce_node_number="<<activenodes<<" "<<cce_node_number<<'\n';
     }
     if(activenodes!=cce_node_number){
       std::cout<<"ERROR: The activenode number of part1 doesn't equal to cce_node_number for part3."<<'\n';
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

     // it is written in DistributePoints and
     // read in DistriPart3zcoords
     mylk0=new LO[li0]; 
     mylk1=new LO[li0];
     mylk2=new LO[li0];      
     DistriPart3zcoords(p1pp3d, test_dir);

// for debugging
     bool debug = false;
     if(debug){
       for(LO i=0;i<li0;i++){
	 LO num=0;
	 LO* recvcount=new LO[p1pp3d.npz];
	 MPI_Allgather(&mylk0[i],1,MPI_INT,recvcount,1,MPI_INT,p1pp3d.comm_z);
	 for(LO h=0;h<p1pp3d.npz;h++) num+=recvcount[h];
	 free(recvcount);
	 if(p1pp3d.mype==2 || p1pp3d.mype==0) std::cout<<num<<"   "<<versurf[li1+i]<<'\n';
	if(i==0){
	  std::cout<<"mypez,lk="<<p1pp3d.mype_z<<" "<<mylk0[0]<<" "<<mylk1[0]<<" "<<mylk2[0]<<'\n';
	} 
      }
    } 
  }
  PERFSTUBS_STOP_STRING(__func__);
}

void Part3Mesh3D::BlockIndexes(const MPI_Comm comm_x,const LO mype_x,const LO npx)
{
  PERFSTUBS_START_STRING(__func__);
  GO* inds = new GO[npx]; 
  blockcount=0;
  for(LO i=0;i<li0;i++)
    blockcount+=(GO)versurf[li1+i];
  inds[mype_x]=blockcount;
  MPI_Datatype mpitype;
  mpitype = getMpiType(GO());
  MPI_Allgather(MPI_IN_PLACE,1,mpitype,inds,1,mpitype,comm_x);
  blockstart=0;
  for(LO i=0;i<mype_x;i++) 
     blockstart+=inds[i];
  blockend=blockstart+blockcount-1;
  delete[] inds;
  PERFSTUBS_STOP_STRING(__func__);
}


void InitzcoordsInCoupler(double* zcoords,LO* versurf,LO nsurf)
{
  PERFSTUBS_START_STRING(__func__);
  double shift=0.1;
  GO num=0;
  for(LO i=0;i<nsurf;i++){
    double delta=2.0*cplPI/(double)(versurf[i]);
    for(LO j=0;j<versurf[i];j++){
      zcoords[num]=(double)(j)*delta+shift;
      num++;
    }
  } 
  PERFSTUBS_STOP_STRING(__func__);
}

//when prepro=ture
void Part3Mesh3D::DistriPart3zcoords(const Part1ParalPar3D &p1pp3d,
    const std::string test_dir)
{
  PERFSTUBS_START_STRING(__func__);
  if(preproc==true){
    if(test_case==TestCase::t0){
      zcoordall = new double[activenodes];
      InitzcoordsInCoupler(zcoordall,versurf,p1pp3d.nx0);
    }else{
      assert(zcoordall); // the number of elements is activenode
    }

    LO numvert=0, numsurf=0;
    for(LO i=0;i<p1pp3d.mype_x;i++){
      for(LO j=xboxinds[i][1];j<xboxinds[i][2]+1;j++){
        numvert+=versurf[numsurf];
	numsurf+=1; 
      } 
    }

    LO index1=p1pp3d.li1;
    LO index2=p1pp3d.li2;
    LO index0=p1pp3d.li0;
    pzcoords = new double*[index0]; 
    zcoordsurf = new double*[index0];
    double* zcoords;
    bool debug=false;
    for(LO i= index1;i<index2+1;i++)
    {
      zcoordsurf[i-index1]=new double[versurf[numsurf]];
      zcoords=zcoordsurf[i-index1];
      for(LO j=0;j<versurf[numsurf];j++){
        zcoords[j]=zcoordall[numvert+j];  
      }
      if(test_case==TestCase::t0){
        assert(!test_dir.empty());
        std::string fname=test_dir+std::to_string(i)+"_zcoords.txt";
        OutputtoFile(zcoords,versurf[numsurf],fname);
      }
      nstart[i] = minloc(zcoords,versurf[numsurf]);
      reshuffleforward(zcoords,nstart[i],versurf[numsurf]);
      if(debug){
	if(p1pp3d.mype==0 && i==index1+1){
	  for(LO h=0;h<versurf[numsurf];h++){
	    std::cout<<"h,pz="<<h<<" "<<zcoordsurf[i-index1][h]<<'\n';
	  }
	}
      }
      DistributePoints(zcoords,index1,i,p1pp3d.pzcoords,p1pp3d);
      pzcoords[i-index1]= new double[mylk0[i-index1]];
      for(LO k=0;k<mylk0[i-index1];k++){
	pzcoords[i-index1][k]= zcoords[mylk1[i-index1]+k];
      } 
      numvert+=versurf[numsurf];
      numsurf+=1;        
    }
  }  
  PERFSTUBS_STOP_STRING(__func__);
}

void Part3Mesh3D::JugeFirstSurfaceMatch(double xp1)
{
  PERFSTUBS_START_STRING(__func__);
  double* tmp = new double[nsurf];
  for(LO i=0;i<nsurf; i++){
    tmp[i]=abs(xcoords[i]-(xp1));
  }
  if(test_case==TestCase::t0){
    shiftx = minloc(tmp,nsurf)-1;
  }else{};
  std::cout<<"cce_first_surface="<<cce_first_surface<<'\n';
  if(shiftx+1!=cce_first_surface){
    std::cout<<"shiftx="<<shiftx<<'\n';
    std::cout<<"ERROR: The first surface of part1 doesn't match cce_first_surface."<<'\n';
    exit(1);
  }
  delete[] tmp;
  PERFSTUBS_STOP_STRING(__func__);
}

LO  minloc(const double* array, const LO n)
{
    PERFSTUBS_START_STRING(__func__);
    double zmin=minimalvalue(array, n);
    LO num=0;
    for(LO i=0;i<n;i++){ 
      num=i;
      if(array[i]==zmin) break;
    }
    PERFSTUBS_STOP_STRING(__func__);
    return num;
 }

//// notice: be carefull with extra_zero case.
// CWS - one of the classes must be read only.... 
void Part3Mesh3D::DistributePoints(const double* exterarr, const LO gstart,LO li, 
                  const double* interarr, const Part1ParalPar3D  &p1pp3d)
{
  PERFSTUBS_START_STRING(__func__);
  if(preproc==true){
    LO nstart;
    double* tmp=new double[versurf[li]];
    for(LO i=0;i<versurf[li];i++){
      tmp[i]=abs(exterarr[i]-interarr[p1pp3d.lk1]);
    }
    nstart=minloc(tmp,versurf[li]);

    //nstart must be in my domain or will duplicate
  //  if(nstart!=0) nstart=nstart-1;
    if(exterarr[nstart]<interarr[p1pp3d.lk1]){
      nstart+=1;
    }
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
  PERFSTUBS_STOP_STRING(__func__);
}

 double minimalvalue(const double* array, const LO n)
{
    PERFSTUBS_START_STRING(__func__);
    double tmp=array[0];
    for(LO i=1;i<n;i++){
      if(tmp>array[i]){
        tmp=array[i];
      }
    }
    PERFSTUBS_STOP_STRING(__func__);
    return tmp;      
}


}
