#include "importpart3mesh.h"
#include "commpart1.h"
#include "testutilities.h"
#include "sendrecv_impl.h"
#include "interpoutil.h"
//#include "macros.h"
#include <cstdlib>
#include <cassert>
#include <cmath>

namespace coupler{

/*|>Initialize XGC's mesh for gene-xgc coupling.*/ 
void Part3Mesh3D::init(const Part1ParalPar3D &p1pp3d,
    const std::string test_dir)
{
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

     /*it is written in DistributePoints and
      * read in DistriPart3zcoords
      */
     mylk0=new LO[li0]; 
     mylk1=new LO[li0];
     mylk2=new LO[li0];      
     DistriPart3zcoords(p1pp3d, test_dir);

     /*for debugging*/
     MPI_Barrier(MPI_COMM_WORLD);
     bool debug = true;
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
}

void Part3Mesh3D::BlockIndexes(const MPI_Comm comm_x,const LO mype_x,const LO npx)
{
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
	if(i==index1+1){
	  for(LO h=0;h<versurf[numsurf];h++){
	    std::cout<<"h,pz="<<h<<" "<<zcoordsurf[i-index1][h]<<'\n';
	  }
	}
      }
      DistributePoints(zcoords, index1, i, p1pp3d.pzcoords, p1pp3d.lk1, p1pp3d.lk2, p1pp3d.nz0);
      pzcoords[i-index1]= new double[mylk0[i-index1]];
      for(LO k=0;k<mylk0[i-index1];k++){
	pzcoords[i-index1][k]= zcoords[mylk1[i-index1]+k];
      } 
      numvert+=versurf[numsurf];
      numsurf+=1;        
    }
  }  
}

void Part3Mesh3D::JugeFirstSurfaceMatch(double xp1)
{
  double* tmp = new double[nsurf];
  for(LO i=0;i<nsurf; i++){
    tmp[i]=xcoords[i]-(xp1);
    if(tmp[i]<=0.0) tmp[i] = -tmp[i];
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
}

LO  minloc(const double* array, const LO n)
{
    double zmin=minimalvalue(array, n);
    LO num=0;
    for(LO i=0;i<n;i++){ 
      num=i;
      if(array[i]==zmin) break;
    }
    return num;
 }

 /*notice: be carefull with extra_zero case.
  *CWS - one of the classes must be read only.... 
  */
void Part3Mesh3D::DistributePoints(const double* exterarr, const LO gstart,LO li, 
                  const double* interarr, const LO lk1, const LO lk2, const LO nz0)
{
  if(preproc==true){
    LO nstart;
    double* tmp=new double[versurf[li]];
    for(LO i=0;i<versurf[li];i++){
      tmp[i]=std::abs(exterarr[i]-interarr[lk1]);
    }
    nstart=minloc(tmp,versurf[li]);

    /*nstart must be in my domain or will duplicate
     */
    if(exterarr[nstart] < interarr[lk1]){
      nstart+=1;
    }
    LO i1=nstart;
    LO i2=nstart;
    double internal_ub;
    if(lk2 == nz0-1){
      internal_ub=cplPI;
    }
    else{
      internal_ub = interarr[lk2 + 1];
    }
    bool inside = true;
    while (inside){
      if (i2>=versurf[li]-1){
        break;
      }
      if (exterarr[i2+1]<internal_ub){
        i2 += 1;
      }
      else{
        inside=false;
      }
    }
    mylk1[li-gstart]=i1;
    mylk2[li-gstart]=i2;
    mylk0[li-gstart]=i2-i1+1;
    if(test_case==TestCase::t0){   
      LO rank;
      MPI_Comm_rank(MPI_COMM_WORLD, &rank);
      std::cout<<"rank="<< rank <<" "<< li-gstart <<'\n';
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

//Initialize XGC's mesh for gem-xgc coupling
void Part3Mesh3D::initXgcGem(const Array1d<LO>* xgccouple,const Array1d<double>* rzcoords)
{
  li0 = p1->li0;
  li1 = p1->li1;
  li2 = p1->li2;
  lj0 = p1->lj0;
  nphi = p1->nphi;
  nwedge = p1->nwedge; 
 
  LO* inttmp;
  inttmp=xgccouple->data();
  totnodes=inttmp[0];
  npsi_surf=inttmp[1];
  cce_first_surface=inttmp[2]; // The number labeling the first active surface 
  cce_last_surface=inttmp[3];
  cce_first_node=inttmp[4];
  cce_last_node=inttmp[5];
  cce_first=inttmp[6];  // The number labeling the first surface; cce_first = 1;
  nsurf=cce_last_surface-cce_first_surface+1;

  if (!p1->mype) fprintf(stderr,"tot:%d, first_node:%d, last_node:%d nsurf: %d, cce_first_surface: %d, cce_last_surface: %d \n",
      totnodes,cce_first_node,cce_last_node, nsurf, cce_first_surface, cce_last_surface);

  versurf = new LO[nsurf];
  nstart  = new LO[nsurf];

  mylk0=new LO[li0]; 
  mylk1=new LO[li0];
  mylk2=new LO[li0];  
  activenodes=0;

  for (LO i=0;i<nsurf;i++){ 
    versurf[i] = inttmp[i + cce_first_surface + 7 - cce_first];  
    activenodes += versurf[i];
  }
  LO totnum = 0;
  for (LO i=0; i<npsi_surf; i++) totnum += inttmp[i + 7];  
  LO begin_node = 0;
  for (LO i = cce_first; i<cce_first_surface; i++) begin_node += inttmp[i + 7 - cce_first];

  if (p1->mype == 0) printf("versurf[0]: %d, versurf[nsurf-1]: %d \n", versurf[0], versurf[nsurf-1]); 
  if (p1->mype == 0) printf("cce_first_node:%d, cce_last_node: %d, activenodes: %d, totnum: %d, begin_node: %d \n", 
     cce_first_node, cce_last_node, activenodes, totnum, begin_node);
  if (!p1->mype == 0) printf("npsi_surf: %d, nsurf: %d totnum-begin_node: %d \n", npsi_surf, 
     nsurf, totnum-begin_node);

  BlockIndexes(p1->comm_x,p1->mype_x,p1->npx); 
  double* realtmp;
  Rcoordall = new double[activenodes];
  Zcoordall = new double[activenodes];
  realtmp = rzcoords->data(); 
  for(LO i=begin_node; i<begin_node+activenodes; i++)  Rcoordall[i - begin_node] = realtmp[i]; 
  for(LO i=totnum+begin_node; i<totnum+begin_node+activenodes; i++){ 
     Zcoordall[i - totnum - begin_node] = realtmp[i];
  }
  if (!p1->mype) printf("Rcoordall[activenodes-1]: %f, Zcoordall[activenodes-1]: %f, %f, %f \n", 
     Rcoordall[activenodes-1], Zcoordall[activenodes-1], Rcoordall[0], Zcoordall[0]);
 
  double eq_axis_r=realtmp[0];
  double eq_axis_z=realtmp[totnodes];
  double* tmp=new double[activenodes];
  for (LO i=0;i<activenodes;i++) 
    tmp[i]=atan2(Zcoordall[i]-eq_axis_z, Rcoordall[i]-eq_axis_r);
    
  theta_xgc = new double*[p1->li0];
  for (LO i=0; i<p1->li0; i++){
    theta_xgc[i] = new double[versurf[p1->li1+i]];
  }

  GO num=blockstart;
  for (LO i=0; i<p1->li0; i++){
    for (LO k=0; k<versurf[p1->li1+i]; k++) theta_xgc[i][k]=tmp[num+k];
    num+=versurf[p1->li1+i];
  }     
/*
    if(p1->mype == 17){
      for( LO m = 0; m<versurf[p1->li1]; m++) printf("m=%d, versurf=%d, tmptheta[i][m]=%f \n", m, versurf[p1->li1],tmptheta[0][m]);
    }
    MPI_Barrier(MPI_COMM_WORLD);
*/
  for (LO i=p1->li1; i<p1->li2+1; i++){  
    nstart[i] = minloc(theta_xgc[i-p1->li1], versurf[i]);
    reshuffleforward(theta_xgc[i-p1->li1], nstart[i], versurf[i]);
    DistributePoints(theta_xgc[i-p1->li1], p1->li1, i, p1->theta, p1->lk1, p1->lk2, p1->nz0);
  } 

  bool debug = true;
  if (debug){
    printf("mype_x: %d, mype: %d, \n", p1->mype_x, p1->mype);

    if (p1->mype == 3){
      for(LO i=0; i<versurf[p1->li1]; i++ )
       printf("nth: %d, theta: %f, versurf[0]: %d \n", i, theta_xgc[0][i], versurf[p1->li1]);
    } 
  }

  debug = true;
  if(debug){
    LO* recvcount=new LO[p1->npz];
    for(LO i=0;i<li0;i++){
/*
      LO num=0;
      MPI_Allgather(&mylk0[i],1,MPI_INT,recvcount,1,MPI_INT,p1->comm_z);
      for(LO h=0;h<p1->npz;h++) num+=recvcount[h];
      if(p1->mype==2 || p1->mype==0) std::cout<<"num, versurf = "<<num<<"   "<<versurf[li1+i]<<'\n';
*/
     if(p1->mype==17){
       std::cout<<"mypez,lk="<<p1->mype_z<<" "<<mylk0[0]<<" "<<mylk1[0]<<" "<<mylk2[0]<<'\n';
     }
   }
   free(recvcount);
 }

  theta_geo=new double*[p1->li0];
  for(LO i=0;i<p1->li0;i++){
    theta_geo[i]=new double[mylk0[i]];
    for(LO k=0;k<mylk0[i];k++) theta_geo[i][k] = theta_xgc[i][mylk1[i]+k];
  }

  
  double* tmpthetaeq = new double[p1->ntheta+5];
  tmpthetaeq[0] = p1->thetagrideq[0]-2.0*p1->dth;
  tmpthetaeq[1] = p1->thetagrideq[0]-p1->dth;
  tmpthetaeq[p1->ntheta+3] = p1->thetagrideq[p1->ntheta]+p1->dth;
  tmpthetaeq[p1->ntheta+4] = p1->thetagrideq[p1->ntheta]+2.0*p1->dth;
  for (LO k=2; k<p1->ntheta+3; k++) tmpthetaeq[k]=p1->thetagrideq[k-2];

  theta_flx = new double*[p1->li0];
  for(LO i=0; i<p1->li0; i++){
     theta_flx[i] = new double[versurf[p1->li1+i]];
  }
 
  /*Here, the continuous boundary condition is used for the 3rd-order Lagrangain interpolaiton; 
   *It's better to replace it with the cubic spline interpolation
   */
  double* tmpflxeq=new double[p1->ntheta+5];

  /*Here: Linear or high order extrapolation may be better to get the boundary value
   */
  // xgc only sends poloida theta. The following is to get the theta flux on each xgc's theta 
  // via interpolating GEM's theta flux mesh. 
  for (LO i=0; i<p1->li0; i++){ 
    tmpflxeq[0] = p1->thflxeq[i][p1->ntheta-2] - 2.0*cplPI; // reason for "-2": thflx[ntheta] = cplPI, thfx[0] = -cplPI
    tmpflxeq[1] = p1->thflxeq[i][p1->ntheta-1] - 2.0*cplPI; // reason for "-1": thflx[ntheta] = cplPI, thfx[0] = -cplPI
    tmpflxeq[p1->ntheta+3] = p1->thflxeq[i][1] + 2.0*cplPI;
    tmpflxeq[p1->ntheta+4] = p1->thflxeq[i][2] + 2.0*cplPI;
    for (LO k=2; k<p1->ntheta+3; k++) tmpflxeq[k] = p1->thflxeq[i][k-2];
    if (theta_xgc[i][0] < tmpthetaeq[1] || theta_xgc[i][p1->li1+i] > tmpthetaeq[p1->ntheta+3]){
       printf("Error: The boundary condition is not right for the Lagrangian interpolation");
       exit(1);
    }   
    Lag3dArray(tmpflxeq,tmpthetaeq,p1->ntheta+5,theta_flx[i], theta_xgc[i],versurf[p1->li1+i]);     
    debug = false;
    if(debug){
      if (p1->mype ==3 && i==0){
	 for (LO j = 0; j < p1->ntheta+5; j++) 
	   fprintf(stderr, "j: %d, tmpflxeq[j]: %f, tmpthetaeq: %f \n", j, tmpflxeq[j], tmpthetaeq[j]);
	 for (LO k = 0; k< versurf[p1->li1]; k++) 
	   fprintf(stderr, "k:%d, theta_flx[0][k]: %f, theta_xgc: %f \n", k, theta_flx[0][k], theta_xgc[0][k]);
      }
    }
  } 
  
  debug = false;
  if(debug){
    if (p1->mype ==3){
       for (LO i = 0; i< versurf[p1->li1]; i++) printf("i:%d, theta_flx[0][i]: %f \n", i, theta_flx[0][i]);
    }
  }

  y_xgc = new double**[p1->li0];
  for(LO i=0;i<p1->li0;i++){
    y_xgc[i] = new double*[p1->nphi];
    for(LO j=0;j<p1->nphi;j++)
      y_xgc[i][j]=new double[mylk0[i]];  
  }
  
  double phi_tmp;
  for(LO i=0;i<p1->li0;i++){
    for(LO j=0;j<p1->nphi;j++){
      phi_tmp=double(j-1)*2.0*cplPI/double(p1->nphi*p1->nwedge);
      for(LO k=0; k<mylk0[i]; k++){
        y_xgc[i][j][k]=remainder(p1->r0/p1->q0*(theta_flx[i][mylk1[i]+k]-phi_tmp),p1->ly);
        if(y_xgc[i][j][k]<0) y_xgc[i][j][k]=p1->ly+y_xgc[i][j][k];
      }
    }
  } 

  theta_pot = new double****[p1->li0];
  theta_ind_pot = new LO****[p1->li0];
  zeta_pot = new double***[p1->li0];
  nodesdist_fl = new double***[p1->li0];
  for(LO i=0;i<p1->li0;i++){
    theta_pot[i] = new double***[p1->lj0];
    theta_ind_pot[i] = new LO***[p1->lj0];
    zeta_pot[i] = new double**[p1->lj0];
    nodesdist_fl[i] = new double**[p1->lj0];
    for(LO j=0;j<p1->lj0;j++){ 
      theta_pot[i][j] = new double**[mylk0[i]];
      theta_ind_pot[i][j] = new LO**[mylk0[i]];
      zeta_pot[i][j] = new double*[mylk0[i]];
      nodesdist_fl[i][j] = new double*[mylk0[i]];
      for(LO k=0;k<mylk0[i];k++){
        theta_pot[i][j][k] = new double*[4];  //For 3rd central Lagrnagian interpolation 
        theta_ind_pot[i][j][k] = new LO*[4]; //For 3rd central Lagrnagian interpolation
        zeta_pot[i][j][k] = new double[5];
        nodesdist_fl[i][j][k] = new double[5];
        for(LO h=0;h<4;h++){
          theta_pot[i][j][k][h] = new double[5];
          theta_ind_pot[i][j][k][h] = new LO[4];
        } 
      }
    }
  }

  theta_ind_gemxgc = new LO**[p1->li0];
  theta_gemxgc = new double**[p1->li0];
  for (LO i=0; i<p1->li0; i++) {
    theta_ind_gemxgc[i] = new LO*[p1->lk0];
    theta_gemxgc[i] = new double*[p1->lk0];
    for (LO k=0; k<p1->lk0; k++) {
      theta_ind_gemxgc[i][k] = new LO[4];
      theta_gemxgc[i][k] = new double[5];
    }
  }

  //for nzb=2;
  double dzeta = 2.0*cplPI/double(p1->nphi*p1->nwedge);
  double lzeta = 2.0*cplPI/double(p1->nwedge);
  double zeta_tmp;
  double y_tmp;
  double q_local;
  double tmptheta;
//  flxxgc flxinter;
  LO j10;

  debug = false;
  for(LO i = 0; i<p1->li0; i++){
    q_local = p1->q_prof[i + p1->li1]; // FIXME 
    for(LO k = 0; k < mylk0[i]; k++){
      y_tmp = double(k)*p1->dy;
      // zeta_tmp corresponding to y given xgc's poloidal theta flux point
      zeta_tmp = remainder(q_local*theta_flx[i][mylk1[i]+k]-p1->q0*y_tmp/p1->r0,2.0*cplPI/double(p1->nwedge));
      if(zeta_tmp < 0) zeta_tmp = zeta_tmp + 2.0*cplPI/double(p1->nwedge);	

      for(LO j = 0; j < p1->nphi; j++){
        zeta_pot[i][j][k][4] = zeta_tmp;
        //FIXME: here is different from GEM
        // zeta_pot stores the 4 zeta points of xgc zeta mesh for the interpolation of zeta_tmp
	j10 = search_zeta(dzeta, lzeta, p1->nphi, zeta_tmp);
        zeta_pot[i][j][k][0] = double(j10-1)*dzeta;
        zeta_pot[i][j][k][1] = double(j10)*dzeta;
        zeta_pot[i][j][k][2] = double(j10+1)*dzeta;
        zeta_pot[i][j][k][3] = double(j10+2)*dzeta;

        debug = false;
	if (debug){
	  if(p1->mype == 3){
	    fprintf(stderr, "theta_flx: %f, zeta_pot: %f, q_local: %f \n", theta_flx[i][mylk1[i] + k], 
	       zeta_pot[i][j][k][0], q_local);
	  }
	}

        // tmpflx is the theta angle  with given theta_flx (the xgc's theta_flx point) and y
        tmptheta = (q_local*theta_flx[i][mylk1[i] + k] - zeta_tmp+zeta_pot[i][j][k][0])/q_local;
      
    //    fprintf(stderr, "tmptheta=: %f \n", tmptheta);        
        search_theta_3rdorder_periodbound(tmptheta, theta_xgc[i], versurf[p1->li1 + i]);
	debug = false;
	if (debug && p1->mype == 1){     
	  printf("flxinter.flxt[i]: %f, %f, %f, %f, %f \n", flxinter.flxt[0], flxinter.flxt[1], flxinter.flxt[2],
	    flxinter.flxt[3], flxinter.flxt[4]);
	  printf("flxinter.flxtind[i]:%d, %d, %d, %d \n", flxinter.flxtind[0], flxinter.flxtind[1], flxinter.flxtind[2], 
	    flxinter.flxtind[3]);
	}

	for (LO h=0; h<5; h++) theta_pot[i][j][k][0][h] = flxinter.flxt[h];
	for (LO h=0; h<4; h++) theta_ind_pot[i][j][k][0][h] = flxinter.flxtind[h];

	tmptheta = (q_local*theta_flx[i][mylk1[i] + k] - zeta_tmp + zeta_pot[i][j][k][1])/q_local;
	search_theta_3rdorder_periodbound(tmptheta,theta_xgc[i],versurf[p1->li1 + i]);
	for (LO h=0; h<5; h++) theta_pot[i][j][k][1][h] = flxinter.flxt[h];
	for (LO h=0; h<4; h++) theta_ind_pot[i][j][k][1][h] = flxinter.flxtind[h];
	debug = false;
	if (debug){     
	  printf("flxinter.flxt[i]: %f, %f, %f, %f, %f \n", flxinter.flxt[0], flxinter.flxt[1], flxinter.flxt[2],
	    flxinter.flxt[3], flxinter.flxt[4]);
	  printf("flxinter.flxtind[i]: %d, %d, %d, %d \n", flxinter.flxtind[0], flxinter.flxtind[1], flxinter.flxtind[2], 
	    flxinter.flxtind[3]);
	}


	tmptheta = (q_local*theta_flx[i][mylk1[i] + k] - zeta_tmp + zeta_pot[i][j][k][2])/q_local;

	search_theta_3rdorder_periodbound(tmptheta, theta_xgc[i], versurf[p1->li1 + i]);
	for (LO h=0; h<5; h++) theta_pot[i][j][k][2][h] = flxinter.flxt[h];
	for (LO h=0; h<4; h++) theta_ind_pot[i][j][k][2][h] = flxinter.flxtind[h];

	tmptheta = (q_local*theta_flx[i][mylk1[i] + k] - zeta_tmp + zeta_pot[i][j][k][3])/q_local;
	search_theta_3rdorder_periodbound(tmptheta, theta_xgc[i], versurf[p1->li1+i]);
	for (LO h=0; h<5; h++) theta_pot[i][j][k][3][h] = flxinter.flxt[h];
	for(LO h=0; h<4; h++)  theta_ind_pot[i][j][k][3][h] = flxinter.flxtind[h];        
        
	nodesdist_fl[i][j][k][0] = 0.0;
	nodesdist_fl[i][j][k][1] = sqrt(pow(theta_pot[i][j][k][1][4] - theta_pot[i][j][k][0][4],2)
	+ pow(zeta_pot[i][j][k][1] - zeta_pot[i][j][k][0],2));

        // this tmptheta is the theta angle of the point where field line with y cross with zeta line of 
        //  xgc labled by theta_flx
        tmptheta = (q_local*theta_flx[i][mylk1[i] + k] - zeta_tmp)/q_local;
        nodesdist_fl[i][j][k][4] = sqrt(pow(tmptheta - theta_pot[i][j][k][0][4],2)
	+ pow(zeta_pot[i][j][k][4] - zeta_pot[i][j][k][0],2));  

        nodesdist_fl[i][j][k][2] = sqrt(pow(theta_pot[i][j][k][2][4] - theta_pot[i][j][k][0][4],2)
	+ pow(zeta_pot[i][j][k][2] - zeta_pot[i][j][k][0],2));    
	nodesdist_fl[i][j][k][3] = sqrt(pow(theta_pot[i][j][k][3][4] - theta_pot[i][j][k][0][4],2)
	+ pow(zeta_pot[i][j][k][3] - zeta_pot[i][j][k][0],2));

        debug = true;
	if (debug && p1->mype == 2 && i==7 && j==15 && k==22){     
/* 
          fprintf(stderr, "fxinter.flxt[i]: %f, %f, %f, %f, %f \n", theta_pot[i][j][k][2][0], theta_pot[i][j][k][2][1],
                 theta_pot[i][j][k][2][2], theta_pot[i][j][k][2][3], theta_pot[i][j][k][2][4], 
                 theta_pot[i][j][k][2][5]);
*/   
//          if(isnan(nodesdist_fl[i][j][k][4])) {
            fprintf(stderr, "nodesdist_fl: %f, %f, %f, %d, %d, %d \n", 
               nodesdist_fl[i][j][k][0], nodesdist_fl[i][j][k][1],  nodesdist_fl[i][j][k][4], i, j, k);
//          }
        }
      }
    }
  }

  MPI_Barrier(MPI_COMM_WORLD);
/*
  fprintf(stderr, "sxz 11 \n");
  for (LO k=0; k<p1->lk0; k++) {
    if (p1->mype == 0) fprintf(stderr, "k: %d, p1->theta: %f \n", k, p1->theta[p1->lk1+k]);  
  }
*/
  InitInterpoThetaPointsXGCtoGEM();

  MPI_Barrier(MPI_COMM_WORLD);

  delete[] tmpthetaeq;
}

void Part3Mesh3D::InitInterpoThetaPointsXGCtoGEM() {  
  for (LO i=0; i<p1->li0; i++ ) {
    for (LO k=0; k<p1->lk0; k++) {
      if (p1->mype == 0) fprintf(stderr, "i: %d, p1->theta: %f \n", i, p1->theta[p1->lk1+k]);  
//      MPI_Barrier(MPI_COMM_WORLD);
      search_theta_3rdorder_periodbound(p1->theta[p1->lk1+k], theta_xgc[i], versurf[p1->li1 + i]); 
      for (LO h=0; h<5; h++) theta_gemxgc[i][k][h] = flxinter.flxt[h];
      for (LO h=0; h<4; h++) theta_ind_gemxgc[i][k][h] = flxinter.flxtind[h];       
    }
  }  
}


void Part3Mesh3D::search_y(LO j1,LO j2,double w1,double w2,const double dy,const double ly,const double tmp)
{
  j1=int(tmp/dy);
  j2=j1+1;
  w2=(tmp-double(j1)*dy)/dy;
  w1=1.0-w2;
}


//fixme: the equivalent routine in GEM looks not right
inline LO Part3Mesh3D::search_zeta(const double dlength,const double length,const LO nlength,double tmp)
{
  LO j10;
  if(tmp<0.0 || tmp>=length){
     tmp=remainder(tmp,length);
     if(tmp<0.0) tmp=tmp+length;
  }
  j10=floor(tmp/dlength);
  return j10;
/*
  if(j10==0){
    j11=-1;
    j20=j10+1;
    j21=j10+2;
  } else if(j10=nlength-1){
    j11=nlength-2;
    j20=0;
    j21=1;
  } else {
    j11=j10-1;
    j20=j10+1;
    j21=j10+2;
  }
  w2=(tmp-real(j1)*dlength)/dlength;
  w1=1.0-w2;
*/
}

 /*This routine is for choosing the points between which that interpolated locates for the 
   third order lagrangian interpolation.*/ 
 void Part3Mesh3D::search_theta_3rdorder_periodbound(
        double tmpflx, const double* flxin, LO num)
{
  bool debug = false;
//  fprintf(stderr, "sxz33 \n");
//  tmpflx = remainder(tmpflx + cplPI,2.0 * cplPI);
//  fprintf(stderr, "tmpflx_1: %f \n", tmpflx);
//  if (tmpflx < 0) tmpflx = tmpflx + 2.0 * cplPI;
//  fprintf(stderr, "tmpflx_2: %f \n", tmpflx);
//  tmpflx = tmpflx - cplPI;
  
//  fprintf(stderr, "tmpflx: %f \n", tmpflx);
  if (tmpflx >= flxin[num - 1]){
    LO k = 0;
    flxinter.flxt[2] = flxin[k] + 2.0*cplPI;
    while(tmpflx > flxinter.flxt[2]){
      k++;
      flxinter.flxt[2] = flxin[k] + 2.0*cplPI;
    }
    flxinter.flxtind[2] = k;
    flxinter.flxtind[3] = k+1;
    if (k >= 3){
      flxinter.flxtind[1] = k-1;
      flxinter.flxtind[0] = k-2;

      flxinter.flxt[0] = flxin[k-2] + 2.0*cplPI;
      flxinter.flxt[1] = flxin[k-1] + 2.0*cplPI;
      flxinter.flxt[2] = flxin[k] + 2.0*cplPI;
      flxinter.flxt[3] = flxin[k+1] + 2.0*cplPI;
    }
    else if (k == 2){
      flxinter.flxtind[1] = 1;
      flxinter.flxtind[0] = num -1;

      flxinter.flxt[0] = flxin[num - 1];   // This assignment is based on flxin has periodic symmetry.  
      flxinter.flxt[1] = flxin[k-1] + 2.0*cplPI;
      flxinter.flxt[2] = flxin[k] + 2.0*cplPI;
      flxinter.flxt[3] = flxin[k+1] + 2.0*cplPI;
    }
    else if(k == 1) {
      flxinter.flxtind[1] = 0;
      flxinter.flxtind[0] = num - 1;

      flxinter.flxt[0] = flxin[num - 1];
      flxinter.flxt[1] = flxin[0];
      flxinter.flxt[2] = flxin[k] + 2.0*cplPI;
      flxinter.flxt[3] = flxin[k+1] + 2.0*cplPI;
    } 
    else {
      flxinter.flxtind[1] = num - 1;
      flxinter.flxtind[0] = num - 2;

      flxinter.flxt[0] = flxin[num - 2];
      flxinter.flxt[1] = flxin[num - 1];
      flxinter.flxt[2] = flxin[k] + 2.0*cplPI;
      flxinter.flxt[3] = flxin[k+1] + 2.0*cplPI; 
    }

   flxinter.flxt[4]=tmpflx;
   if(debug){
     printf("flxt: %f, %f, %f, %f \n", flxinter.flxt[0], flxinter.flxt[1], flxinter.flxt[2], flxinter.flxt[3]);
   }
   if (tmpflx > flxinter.flxt[2]){
      printf("Error: The boundary condition is not right at upper end first,tmpflx: %f, flxt[2]: %f \n", 
	tmpflx, flxinter.flxt[2]);
      exit(1);
   }

  }
  else if (tmpflx <= flxin[0]){
    LO k = num - 2;
    flxinter.flxt[2] = flxin[k] - 2.0*cplPI;
    while (tmpflx < flxinter.flxt[2]){
      k--;
      flxinter.flxt[2] = flxin[k] - 2.0*cplPI;
    } 
    if (k < num - 2){
      flxinter.flxtind[2] = k;
      flxinter.flxtind[3] = k+1; 
      flxinter.flxtind[1] = k-1;
      flxinter.flxtind[0] = k-2;

      flxinter.flxt[0] = flxin[k-2] - 2.0*cplPI;
      flxinter.flxt[1] = flxin[k-1] - 2.0*cplPI;
      flxinter.flxt[2] = flxin[k] - 2.0*cplPI;
      flxinter.flxt[3] = flxin[k+1] - 2.0*cplPI;
    } 
    else{
      flxinter.flxtind[2] = 0;
      flxinter.flxtind[3] = 1; 
      flxinter.flxtind[1] = num - 2;
      flxinter.flxtind[0] = num - 3;

      flxinter.flxt[0] = flxin[num - 3] - 2.0*cplPI;
      flxinter.flxt[1] = flxin[num - 2] - 2.0*cplPI;
      flxinter.flxt[2] = flxin[0];
      flxinter.flxt[3] = flxin[1]; 
    } 

    flxinter.flxt[4]=tmpflx;
    
    if(debug){
      printf("flxt: %f, %f, %f, %f \n", flxinter.flxt[0], flxinter.flxt[1], flxinter.flxt[2], flxinter.flxt[3]);
    }
    if (tmpflx<=flxinter.flxt[1]){
      printf("Error: The boundary condition is not right at lower end first,tmpflx: %f, flxt[1]: %f \n", 
	tmpflx, flxinter.flxt[1]);
      exit(1);
    }
  }
  else if (tmpflx<=flxin[1] && tmpflx>flxin[0]){
    flxinter.flxt[4]=tmpflx;
    flxinter.flxt[0]=flxin[num-2]-2.0*cplPI;
    flxinter.flxt[1]=flxin[0];
    flxinter.flxt[2]=flxin[1];
    flxinter.flxt[3]=flxin[2];
    
    flxinter.flxtind[0]=num-2;
    flxinter.flxtind[1]=0;
    flxinter.flxtind[2]=1;
    flxinter.flxtind[3]=2;
  }
    

  else if (tmpflx>=flxin[num-2] && tmpflx<flxin[num-1]){
    flxinter.flxt[4]=tmpflx;
    flxinter.flxt[0]=flxin[num-3];
    flxinter.flxt[1]=flxin[num-2];
    flxinter.flxt[2]=flxin[num-1];
    flxinter.flxt[3]=flxin[1]+2.0*cplPI;

    flxinter.flxtind[0]=num-3;
    flxinter.flxtind[1]=num-2;
    flxinter.flxtind[2]=num-1;
    flxinter.flxtind[3]=1;
  }

  else {
    LO k=0;
    while(tmpflx >= flxin[k] && tmpflx != flxin[1] && tmpflx != flxin[num-2]) k+=1;
    if(tmpflx >= flxin[k]) k++;
    if (debug) {
      for (LO i = 0; i< versurf[0]; i++) printf("i:%d, flxin[i]: %f \n", i, flxin[i]);
      assert(num-1 >= k+2);
      fprintf(stderr, "num-1: %d, k-2: %d, k-1: %d, k+1: %d \n", num-1,k-2,k-1,k+1);
      fprintf(stderr, "mype: %d, tmpflx: %f \n", p1->mype, tmpflx);
      fprintf(stderr, "%f, %f, %f, %f, %f \n", tmpflx, flxin[k-2], flxin[k-1], flxin[k], flxin[k+1]);
    }

    flxinter.flxt[4]=tmpflx;
    flxinter.flxt[0]=flxin[k-2];
    flxinter.flxt[1]=flxin[k-1];
    flxinter.flxt[2]=flxin[k];
    flxinter.flxt[3]=flxin[k+1];

    flxinter.flxtind[0]=k-2;
    flxinter.flxtind[1]=k-1;
    flxinter.flxtind[2]=k;
    flxinter.flxtind[3]=k+1; 
    debug = false;  
    if (debug) { 
      fprintf(stderr, "flxinter.flxtind[i]: %d, %d, %d, %d \n", flxinter.flxtind[0], flxinter.flxtind[1], flxinter.flxtind[2], flxinter.flxtind[3]);
    }
  } 

}

} /*importpart3mesh.cc*/ 
