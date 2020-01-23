# default
es: $(OBJ) $(OBJ_ES) 
	$(CMP) -o xgc-es $(OBJ_ES) $(OBJ) $(LIB) ${PETSC_POST_LINK_OPTS} ${PETSC_KSP_LIB} -mkl

em: em_build $(OBJ) $(OBJ_EM) 
	$(CMP) -o xgc-em $(OBJ_EM) $(OBJ) $(LIB) ${PETSC_POST_LINK_OPTS} ${PETSC_KSP_LIB}

es-gpu: gpu_build $(OBJ) $(OBJ_GPU) $(OBJ_ES) 
	$(CMP) -o xgc-es-gpu $(GPU_LINK_OPT) $(OBJ) $(OBJ_GPU) $(OBJ_ES) $(LIB) $(PETSC_KSP_LIB) 

em-gpu: em_build gpu_build $(OBJ) $(OBJ_GPU) $(OBJ_EM) 
	$(CMP) -o xgc-em-gpu $(GPU_LINK_OPT) $(OBJ_EM) $(OBJ_GPU) $(OBJ) $(LIB) 

em_build:
	$(eval FC_FLAGS += -DXGC1_EM)

gpu_build:
	$(eval FC_FLAGS += -DUSE_GPU=1)

$(OBJ): module.F90 search.F90 pol_decomp.F90 f0module.F90
$(OBJ_EM): module.F90 search.F90 pol_decomp.F90 f0module.F90
$(OBJ_ES): module.F90 search.F90 pol_decomp.F90 f0module.F90

clean:
	rm -f *.o *.mod *.lst

.F90.o :
#	-@echo ${FC_FLAGS}
	$(CMP) ${FC_FLAGS} ${FCPPFLAGS} ${XGC_FLAGS} ${XGC_INCLUDE} ${PETSC_INCLUDE_OPTS} -c $<

#
# GPU
#

push_mod_gpu.o:  $(GPU_SRC) push_mod_gpu.F90
	$(CMP) $(XGC_FLAGS) ${GPU_FFLAGS} -Mpreprocess -E push_mod_gpu.F90 > push_mod_gpu.f
	$(CMP) $(XGC_FLAGS) $(GPU_FFLAGS) ${XGC_INCLUDE} -Mfree -c push_mod_gpu.f
