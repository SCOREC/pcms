
/*
 * -- SuperLU routine (version 3.0) --
 * Univ. of California Berkeley, Xerox Palo Alto Research Center,
 * and Lawrence Berkeley National Lab.
 * October 15, 2003
 *
 *
 * -- Modified for XGC0 in Fortran90 --
 *  by E.S. Yoon
 * "iopt" has been removed and unified in that we don't have to distinguish
 * the factorization procedure as individual.
 * June 05, 2011
 */

#include "slu_ddefs.h"

#define HANDLE_SIZE  8
/* kind of integer to hold a pointer.  Use int.
   This might need to be changed on 64-bit systems. */
typedef long long fptr;  /* 64-bit by default */

typedef struct {
    SuperMatrix *L;
    SuperMatrix *U;
    int *perm_c;
    int *perm_r;
} factors_t;

void
c_fortran_dgssv_( int *n, int *nnz, int *nrhs, 
                 double *values, int *rowind, int *colptr,
                 double *b, int *ldb,
		 fptr *f_factors, /* a handle containing the address
				     pointing to the factored matrices */
		 int *info, int *info2)

{
    SuperMatrix A, AC, B;
    SuperMatrix *L, *U;
    int *perm_r; /* row permutations from partial pivoting */
    int *perm_c; /* column permutation vector */
    int *etree;  /* column elimination tree */
    SCformat *Lstore;
    NCformat *Ustore;
    int      i, panel_size, permc_spec, relax;
    trans_t  trans;
    mem_usage_t   mem_usage;
    superlu_options_t options;
    SuperLUStat_t stat;

    trans = NOTRANS;

    /* LU decomposition */

    /* Set the default input options. */
    set_default_options(&options);

	/* Initialize the statistics variables. */
	StatInit(&stat);

	/* Adjust to 0-based indexing */
	for (i = 0; i < *nnz; ++i) --rowind[i];
	for (i = 0; i <= *n; ++i) --colptr[i];

	dCreate_CompCol_Matrix(&A, *n, *n, *nnz, values, rowind, colptr,
			       SLU_NC, SLU_D, SLU_GE);
	L = (SuperMatrix *) SUPERLU_MALLOC( sizeof(SuperMatrix) );
	U = (SuperMatrix *) SUPERLU_MALLOC( sizeof(SuperMatrix) );
	if ( !(perm_r = intMalloc(*n)) ) ABORT("Malloc fails for perm_r[].");
	if ( !(perm_c = intMalloc(*n)) ) ABORT("Malloc fails for perm_c[].");
	if ( !(etree = intMalloc(*n)) ) ABORT("Malloc fails for etree[].");

	/*
	 * Get column permutation vector perm_c[], according to permc_spec:
	 *   permc_spec = 0: natural ordering 
	 *   permc_spec = 1: minimum degree on structure of A'*A
	 *   permc_spec = 2: minimum degree on structure of A'+A
	 *   permc_spec = 3: approximate minimum degree for unsymmetric matrices
	 */    	
	permc_spec = options.ColPerm;
	get_perm_c(permc_spec, &A, perm_c);
	
	sp_preorder(&options, &A, perm_c, etree, &AC);

	panel_size = sp_ienv(1);
	relax = sp_ienv(2);

	dgstrf(&options, &AC, relax, panel_size, etree,
                NULL, 0, perm_c, perm_r, L, U, &stat, info);

	if ( *info == 0 ) {
	    Lstore = (SCformat *) L->Store;
	    Ustore = (NCformat *) U->Store;
#if defined(VPIC_DEBUG_SPARSE)	
	    printf("No of nonzeros in factor L = %d\n", Lstore->nnz);
	    printf("No of nonzeros in factor U = %d\n", Ustore->nnz);
	    printf("No of nonzeros in L+U = %d\n", Lstore->nnz + Ustore->nnz);
#endif
	    dQuerySpace(L, U, &mem_usage);
#if defined(VPIC_DEBUG_SPARSE)		    
	    printf("L\\U MB %.3f\ttotal MB needed %.3f\n",
		   mem_usage.for_lu/1e6, mem_usage.total_needed/1e6);
#endif		   
	} else {
	    printf("dgstrf() error returns INFO= %d\n", *info);
	    if ( *info <= *n ) { /* factorization completes */
		dQuerySpace(L, U, &mem_usage);
		printf("L\\U MB %.3f\ttotal MB needed %.3f\n",
		       mem_usage.for_lu/1e6, mem_usage.total_needed/1e6);
	    }
    }	
	
	/* Restore to 1-based indexing */
	for (i = 0; i < *nnz; ++i) ++rowind[i];
	for (i = 0; i <= *n; ++i) ++colptr[i];

	/* Free un-wanted storage */
	SUPERLU_FREE(etree);
	Destroy_SuperMatrix_Store(&A);
	Destroy_CompCol_Permuted(&AC);
	StatFree(&stat);

    /* Triangular solve */
	/* Initialize the statistics variables. */
	StatInit(&stat);

	dCreate_Dense_Matrix(&B, *n, *nrhs, b, *ldb, SLU_DN, SLU_D, SLU_GE);

    /* Solve the system A*X=B, overwriting B with X. */
    dgstrs (trans, L, U, perm_c, perm_r, &B, &stat, info2);

    Destroy_SuperMatrix_Store(&B);
    StatFree(&stat);

    SUPERLU_FREE (perm_r);
    SUPERLU_FREE (perm_c);
    Destroy_SuperNode_Matrix(L);
    Destroy_CompCol_Matrix(U);
    SUPERLU_FREE (L);
    SUPERLU_FREE (U);
}


