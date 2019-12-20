/* cbicub.cu - bicubic spline evaluation using textures
*
*
*  2D cuda arrays and texture memory have the advantage of block
*  linear addressing that give 2D data locality, potentially reducing
*  the number of fetches if we traverse the grid correctly.
*
*  I'm going to see if that speeds up the critical bicub_spline
*  interpolation routine.
*
*
*/

extern "C" {
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
// Includes CUDA
#include <cuda_runtime.h>
#include <cuda.h>

struct bicub_obj {
	int nr, nz, ncoeff;
	double rmin, dr_inv, zmin, dz_inv;
    cudaArray *rc_cub;
	cudaArray *zc_cub;
	cudaArray *acoeff_cub;
};

texture<int2, cudaTextureType1D, cudaReadModeElementType> tex_rc;
texture<int2, cudaTextureType1D, cudaReadModeElementType> tex_zc;
texture<int2, cudaTextureType2DLayered, cudaReadModeElementType> tex_acoeff;

// FIXME :: Hardcoding ncoeff == 4!!
// I really hate static persistant objects, but I can't really think
// of an easier way to do this right now.

static struct bicub_obj *host_psi_bco;
__device__ __constant__ struct bicub_obj d_psi_bco;


#define gpuErrchk(ans) { gpuAssert((ans), __FILE__, __LINE__); }
inline void gpuAssert(cudaError_t code, const char *file, int line, bool abort=true)
{
   if (code != cudaSuccess) 
   {
      fprintf(stderr,"GPUassert: %s %s %d\n", cudaGetErrorString(code), file, line);
      if (abort) exit(code);
   }
}

// Round 1 - allocate spline memory in C, and use it as normal.

#define FINDEX(bco, c1, c2, ri, zi) (4 * ( 4 * ( (bco).nr * (zi) + (ri)) + (c2)) + (c1))
#define CINDEX(bco, c1, c2, ri, zi) ((bco).nr * ( (bco).nz * ( 4 * (c2) + (c1)) + (zi)) + (ri))
#define FFLAT(bco, c, ri, zi) (16 * ( (bco).nr * (zi) + (ri)) + (c) )

#define COEFFBLOCK(c1, c2) (4 * (c2) + (c1))
// init_cuda_bicub
//
// Allocate memory for the bicubic spline coefficients, and copy them.
void
init_cuda_bicub(double *F_rc_cub, double *F_zc_cub, double *F_acoeff,
	int ncoeff, int nr, int nz,
	double rmin, double zmin, double dr_inv, double dz_inv)
{

	if (host_psi_bco) return;
	host_psi_bco = (struct bicub_obj *) calloc(1, sizeof(bicub_obj));

	host_psi_bco->nr = nr;
	host_psi_bco->nz = nz;
	host_psi_bco->ncoeff = ncoeff;	
	host_psi_bco->rmin	= rmin;
	host_psi_bco->zmin = zmin;
	host_psi_bco->dr_inv = dr_inv;
	host_psi_bco->dz_inv = dz_inv;

	// Create a channel with 2 int, which will store two halves of the double
    cudaChannelFormatDesc channelDesc = cudaCreateChannelDesc(32, 32, 0, 0, cudaChannelFormatKindSigned);

    // ================================
    // Setup and bind the 1D rc/zc arrays
    // ================================

	// Simple standard cuda arrays right now
	gpuErrchk(cudaMallocArray(&(host_psi_bco->rc_cub), &channelDesc, nr));
	gpuErrchk(cudaMallocArray(&(host_psi_bco->zc_cub), &channelDesc, nz));		


	// Copy straight out of Fortran without any reshaping
	gpuErrchk(cudaMemcpyToArray(host_psi_bco->rc_cub, 0, 0, (void *)F_rc_cub,
		sizeof(double)*nr, cudaMemcpyHostToDevice));

	gpuErrchk(cudaMemcpyToArray(host_psi_bco->zc_cub, 0, 0, (void *)F_zc_cub,
		sizeof(double)*nz, cudaMemcpyHostToDevice));

	// Bind the textures to the arrays, with zeros if we go outside, and accessing
	// through un-normalized coordinated with no interpolation.
    tex_rc.addressMode[0] = cudaAddressModeClamp;
    tex_rc.filterMode = cudaFilterModePoint;
    tex_rc.normalized = false;    // access with normalized texture coordinates
    gpuErrchk(cudaBindTextureToArray(tex_rc, host_psi_bco->rc_cub,channelDesc));

    tex_zc.addressMode[0] = cudaAddressModeClamp;
    tex_zc.filterMode = cudaFilterModePoint;
    tex_zc.normalized = false;    // access with normalized texture coordinates
    gpuErrchk(cudaBindTextureToArray(tex_zc, host_psi_bco->zc_cub,channelDesc));

    // ================================
    // Setup and bind the 2D layered coefficient array
    // ================================

    gpuErrchk(cudaMalloc3DArray(&(host_psi_bco->acoeff_cub), &channelDesc, 
    	make_cudaExtent(nr, nz, ncoeff*ncoeff), cudaArrayLayered));

    // I think we need to reshape the fortran array..
    double re_coeff[nr*nz*ncoeff*ncoeff];
    for (int ir=0; ir < nr; ir++) {
    	for (int iz=0; iz < nz; iz++) {
    		for (int c2=0; c2 < ncoeff; c2++) {
    			for (int c1=0; c1 < ncoeff; c1++) {
    				re_coeff[CINDEX(*host_psi_bco,c1,c2,ir,iz)] = F_acoeff[FINDEX(*host_psi_bco,c1,c2,ir,iz)];
    			}
    		}
    	}
    }


    cudaMemcpy3DParms myparms = {0};
    myparms.srcPos = make_cudaPos(0,0,0);
    myparms.dstPos = make_cudaPos(0,0,0);
    myparms.srcPtr = make_cudaPitchedPtr(re_coeff, nr * sizeof(double), nr, nz);
    myparms.dstArray = host_psi_bco->acoeff_cub;
    myparms.extent = make_cudaExtent(nr, nz, ncoeff*ncoeff);
    myparms.kind = cudaMemcpyHostToDevice;
    gpuErrchk(cudaMemcpy3D(&myparms));


    // Bind the acoeff texture 
    tex_acoeff.addressMode[0] = cudaAddressModeClamp;
    tex_acoeff.addressMode[1] = cudaAddressModeClamp;
    tex_acoeff.filterMode = cudaFilterModePoint;
    tex_acoeff.normalized = false;
    gpuErrchk(cudaBindTextureToArray(tex_acoeff,host_psi_bco->acoeff_cub,channelDesc));


	gpuErrchk(cudaMemcpyToSymbol(d_psi_bco, host_psi_bco, sizeof(struct bicub_obj)));

}

void
destroy_cuda_bicub(void)
{
	gpuErrchk(cudaUnbindTexture(tex_rc));
	gpuErrchk(cudaUnbindTexture(tex_zc));	
	gpuErrchk(cudaUnbindTexture(tex_acoeff));		
	gpuErrchk(cudaFreeArray(host_psi_bco->rc_cub));
	gpuErrchk(cudaFreeArray(host_psi_bco->zc_cub));	
	gpuErrchk(cudaFreeArray(host_psi_bco->acoeff_cub));		

	free(host_psi_bco);
	host_psi_bco = NULL;
}

// The 'get' uses a floating point index, like the cuda tex** calls
static __inline__ __device__ double get1D_double(texture<int2, cudaTextureType1D> t, float i)
{
int2 v = tex1D(t,i);
return __hiloint2double(v.y, v.x);
}

// The 'get' uses a floating point index, like the cuda tex** calls
static __inline__ __device__ double get2DLayer_double(texture<int2, cudaTextureType2DLayered> t, float i0, float i1, int lay)
{
int2 v = tex2DLayered(t,i0,i1,lay);
return __hiloint2double(v.y, v.x);
}

// The 'fetch' uses an integer index, like the cuda texfetch** calls
static __inline__ __device__ double fetch1D_double(texture<int2, cudaTextureType1D> t, int i)
{
int2 v = tex1Dfetch(t,i);
return __hiloint2double(v.y, v.x);
}

__device__ void
fetch_cuda_bicub(double *F_rc_cub, double* F_zc_cub, double *F_acoeff, int ir, int iz)
{

	*F_rc_cub = get1D_double(tex_rc,(float)(ir-1));
	*F_zc_cub = get1D_double(tex_zc,(float)(iz-1));
	#pragma unroll
	for (int c = 0; c < 16; c++) {
		F_acoeff[c] = get2DLayer_double(tex_acoeff, (float)(ir-1), (float)(iz-1), c);
	}
}

//#define A(c1,c2) d_psi_bco.acoeff_cub[FINDEX(d_psi_bco,(c1),(c2),ix,iy)]
//#define A(c1,c2) fetch1D_double(tex_acoeff,FINDEX(d_psi_bco,(c1),(c2),ix,iy))
#define A(c1,c2) get2DLayer_double(tex_acoeff, f_ix, f_iy, COEFFBLOCK(c1,c2))
__device__ void
eval_0_cuda(double x, double y, double *f00_)
{
	float f_ix = (x - d_psi_bco.rmin) * d_psi_bco.dr_inv;
	float f_iy = (y - d_psi_bco.zmin) * d_psi_bco.dz_inv;	

	double dx = x - get1D_double(tex_rc, f_ix);
	double dy = y - get1D_double(tex_zc, f_iy);;

	double fy_i;

	double f00;

	f00 = 0;

	fy_i = ((A(0,3)*dy + A(0,2))*dy + A(0,1))*dy + A(0,0);
	f00 = f00 + fy_i;

	fy_i = ((A(1,3)*dy + A(1,2))*dy + A(1,1))*dy + A(1,0);
	f00 = f00 + dx*fy_i;

	fy_i = ((A(2,3)*dy + A(2,2))*dy + A(2,1))*dy + A(2,0);
	f00 = f00 + (dx*dx)*fy_i;

	fy_i = ((A(3,3)*dy + A(3,2))*dy + A(3,1))*dy + A(3,0);
	f00 = f00 + (dx*dx)*dx*fy_i;

    *f00_ = f00;
}

__device__ void
eval_1_cuda(double x, double y, double *f00_, double *f10_, double *f01_)
{
	float f_ix = (x - d_psi_bco.rmin) * d_psi_bco.dr_inv;
	float f_iy = (y - d_psi_bco.zmin) * d_psi_bco.dz_inv;	

	double dx = x - get1D_double(tex_rc, f_ix);
	double dy = y - get1D_double(tex_zc, f_iy);;

	double dfx_i, fx_i;

	double f00, f10, f01;

	f00 = 0;
	f01 = 0;

	fx_i = ((A(3,0)*dx + A(2,0))*dx + A(1,0))*dx + A(0,0);
	f00 = f00 + fx_i;

	fx_i = ((A(3,1)*dx + A(2,1))*dx + A(1,1))*dx + A(0,1);
	f00 = f00 + dy*fx_i;
	f01 = f01 + fx_i;

	fx_i = ((A(3,2)*dx + A(2,2))*dx + A(1,2))*dx + A(0,2);
	f00 = f00 + (dy*dy)*fx_i;
	f01 = f01 + 2.0*dy*fx_i;

	fx_i = ((A(3,3)*dx + A(2,3))*dx + A(1,3))*dx + A(0,3);
	f00 = f00 + dy*((dy*dy)*fx_i);
	f01 = f01 + 3.0*((dy*dy)*fx_i);
	
	f10 = 0;

	dfx_i = (A(3,0)*3.0*dx + A(2,0)*2.0)*dx + A(1,0);
	f10 = f10 + dfx_i;

	dfx_i = (A(3,1)*3.0*dx + A(2,1)*2.0)*dx + A(1,1);
	f10 = f10 + dy*dfx_i;

	dfx_i = (A(3,2)*3.0*dx + A(2,2)*2.0)*dx + A(1,2);
	f10 = f10 + (dy*dy)*dfx_i;

	dfx_i = (A(3,3)*3.0*dx + A(2,3)*2.0)*dx + A(1,3);
	f10 = f10 + (dy*dy)*dy*dfx_i;

	*f00_ = f00;
	*f10_ = f10;
	*f01_ = f01;
}

__device__ void
eval_2_cuda(double x, double y, double *f00_,
	double *f10_, double *f01_, 
	double *f11_, double *f20_, double *f02_)
{
	/*
    ! ----------------------------------
    ! evaluate bicubic polynomial f(x,y)
    ! and high order derivatives
    !
    ! note (xc,yc) is offset or center of box
    !
    ! f00 = f(x,y)
    ! f10 = df/dx
    ! f01 = df/dy
    ! f11 = df^2/dx/dy
    ! f20 = df^2/dx/dx
    ! f02 = df^2/dy/dy
    ! ----------------------------------
	*/
	float f_ix = (x - d_psi_bco.rmin) * d_psi_bco.dr_inv;
	float f_iy = (y - d_psi_bco.zmin) * d_psi_bco.dz_inv;	

	double dx = x - get1D_double(tex_rc, f_ix);
	double dy = y - get1D_double(tex_zc, f_iy);;

	double fx_i, dfx_i, dfx2_i;

	double f00, f10, f01, f11,f20,f02;

	f00 = 0;
	f01 = 0;
	f02 = 0;
 
	fx_i = ((A(3,0)*dx + A(2,0))*dx + A(1,0))*dx + A(0,0);
	f00 = f00 + fx_i;

	fx_i = ((A(3,1)*dx + A(2,1))*dx + A(1,1))*dx + A(0,1);
	f00 = f00 + dy*fx_i;
	f01 = f01 +    fx_i;

	fx_i = ((A(3,2)*dx + A(2,2))*dx + A(1,2))*dx + A(0,2);
	f00 = f00 + dy*(dy*fx_i);
	f01 = f01 + 2.0*(dy*fx_i);
	f02 = f02 + 2.0*fx_i;

	fx_i = ((A(3,3)*dx + A(2,3))*dx + A(1,3))*dx + A(0,3);
	f00 = f00 + dy*(dy*(dy*fx_i));
	f01 = f01 + 3.0*(dy*(dy*fx_i));
	f02 = f02 + 6.0*(dy*fx_i);

	f10 = 0;
	f11 = 0;

	dfx_i = (A(3,0)*3.0*dx + A(2,0)*2.0)*dx + A(1,0);
	f10 = f10 + dfx_i;

	dfx_i = (A(3,1)*3.0*dx + A(2,1)*2.0)*dx + A(1,1);
	f10 = f10 + dy*dfx_i;
	f11 = f11 +    dfx_i;

	dfx_i = (A(3,2)*3.0*dx + A(2,2)*2.0)*dx + A(1,2);
	f10 = f10 + dy*(dy*dfx_i);
	f11 = f11 + 2.0*(dy*dfx_i);

	dfx_i = (A(3,3)*3.0*dx + A(2,3)*2.0)*dx + A(1,3);
	f10 = f10 + dy*(dy*dy*dfx_i);
	f11 = f11 + 3.0*(dy*dy*dfx_i);

	f20 = 0;

	dfx2_i = (3.*2.*A(3,0)*dx + 2.*1.*A(2,0));
	f20 = f20 + dfx2_i;

	dfx2_i = (3.*2.*A(3,1)*dx + 2.*1.*A(2,1));
	f20 = f20 + dy*dfx2_i;

	dfx2_i = (3.*2.*A(3,2)*dx + 2.*1.*A(2,2));
	f20 = f20 + (dy*dy)*dfx2_i;

	dfx2_i = (3.*2.*A(3,3)*dx + 2.*1.*A(2,3));
	f20 = f20 + (dy*dy)*dy*dfx2_i;

	*f00_ = f00;
	*f10_ = f10;
	*f01_ = f01;
	*f11_ = f11;
	*f20_ = f20;
	*f02_ = f02;
}

};