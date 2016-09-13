# ===============================================================
# 			UNIVERSIDADE FEDERAL DE PERNANMBUCO
#           DEPARTAMENTO DE ENGENHARIA MECANICA (DEMEC-CTG-UFPE)
#			NUCLEO DE TECNOLOGIA (NT-CAA-UFPE)
#
# Authors: 	Rogerio Soares(rogerio.soaress@ufpe.br)
#           Darlan Carvalho
#			Saulo
#			Guilherme Caminha
#			Paulo Lyra
# Created:	2008-2015
# ===============================================================

# FLAGS LEGENDS:
# 	FVPO                   : For Visualization Purposes Only
# 	-g:                    : Debugging compilation (Use it to debug program by DDD)
# 	__ADAPTATION_DEBUG__   : Debug Adaptation functions only
#	__ERROR_ANALYSIS_DEBUG__: Debug error analysis functions only
#	_SEEKFORBUGS_           : Debuggin in general (whole code)
#    NOADAPTATION			: Compile code without mesh adaptation feature
#

CXXFLAGS=-g -Wall -DNOADAPTATION -DTRACKING_PROGRAM_STEPS -Wparentheses -D_SEEKFORBUGS_ 
#CXXFLAGS=-g -DNOADAPTATION -DPARALLEL -DTRACKING_PROGRAM_STEPS
#CXXFLAGS=-DPARALLEL -g -DNOADAPTATION -Wunused-local-typedefs -D_SEEKFORBUGS_ -DTRACKING_PROGRAM_STEPS
#CXXFLAGS=-DPARALLEL -g -Wall -Wunused -D__ADAPTATION_DEBUG__ -D__ERROR_ANALYSIS_DEBUG__ -D_SEEKFORBUGS_ -DTRACKING_PROGRAM_STEPS -DFVPO
# scientific_application folder is where libraries are installed 
APP_DIR=$(HOME)/applications
PROJ_DIR=$(HOME)/interpolation/padmec-amr
GMSH_DIR=$(APP_DIR)/gmshGMSH

# compilers
CXX=mpicxx

# paths for headers
INCLUDES=-I$(PROJ_DIR)/include \
	-I$(PROJ_DIR)/src/adaptation -I$(PROJ_DIR)/src/adaptation/adaptive-remeshing -I$(PROJ_DIR)/src/adaptation/h-refinement \
	-I$(PROJ_DIR)/src/adaptation/rh-refinement -I$(APP_DIR)/MAdLib_no_gmsh/include \
	-I$(PROJ_DIR)/src/elliptic -I$(PROJ_DIR)/src/elliptic/MEBFV -I$(PROJ_DIR)/src/hyperbolic \
	-I$(PROJ_DIR)/src/error-analysis -I$(PROJ_DIR)/src/interpolation -I$(PROJ_DIR)/src/pre-processor/EBFV1 -I$(PROJ_DIR)/src/pre-processor/EBFV1_modified \
	-I$(PROJ_DIR)/src/SIMULATION_core \
	-I$(PROJ_DIR)/src/tools -I$(PROJ_DIR)/src/tools/GeomData \
	-I$(GMSH_DIR)/build/Common -I$(GMSH_DIR)/Common -I$(GMSH_DIR)/Geo -I$(GMSH_DIR)/Mesh -I$(GMSH_DIR)/Numeric -I$(GMSH_DIR)/Parser -I$(GMSH_DIR)/Plugin -I$(GMSH_DIR)/Post \
 	-I$(APP_DIR)/FMDB-2011/include -I$(PETSC_DIR)/include -I$(PETSC_DIR)/arch-linux2-c-opt/include 

LIBS1=-L$(APP_DIR)/FMDB-2011/lib -lFMDB-O

ifeq (,$(findstring NOADAPTATION,$(CXXFLAGS) $(ROCKPROP)))
LIBS2=-L$(APP_DIR)/gmshGMSH/lib -lGmsh
endif

LIBS=$(LIBS1) $(LIBS2)

# path for where the objects must be put
# ======================================
OBJ_DIR=$(PROJ_DIR)/objs

# paths for each program source part
# ==================================
SRC_DIR1=$(PROJ_DIR)/src
SRC_DIR2=$(PROJ_DIR)/src/elliptic/EBFV
SRC_DIR21=$(PROJ_DIR)/src/elliptic/MEBFV
SRC_DIR4=$(PROJ_DIR)/src/hyperbolic
SRC_DIR5=$(PROJ_DIR)/src/pre-processor/EBFV1
SRC_DIR51=$(PROJ_DIR)/src/pre-processor/EBFV1_modified
SRC_DIR6=$(PROJ_DIR)/src/SIMULATION_core
SRC_DIR7=$(PROJ_DIR)/src/tools
SRC_DIR71=$(PROJ_DIR)/src/tools/GeomData
SRC_DIR8=$(PROJ_DIR)/src/error-analysis
SRC_DIR9=$(PROJ_DIR)/src/adaptation
SRC_DIR10=$(PROJ_DIR)/src/interpolation
SRC_DIR11=$(PROJ_DIR)/src/adaptation/adaptive-remeshing
SRC_DIR12=$(PROJ_DIR)/src/adaptation/h-refinement

OBJS_MAIN=$(OBJ_DIR)/main.o $(OBJ_DIR)/SIMULATION_core.o $(OBJ_DIR)/SIMULATION_core__solvers.o $(OBJ_DIR)/SIMULATION_adaptation.o \
          $(OBJ_DIR)/CalculateElementsError.o $(OBJ_DIR)/CalculateGlobalError.o $(OBJ_DIR)/ErrorAnalysis.o $(OBJ_DIR)/ErrorAnalysisAuxiliar.o \
          $(OBJ_DIR)/CalculateSmoothGradientNorm.o $(OBJ_DIR)/Calculate_heights.o $(OBJ_DIR)/IhR_main.o $(OBJ_DIR)/interpolation.o \
	      $(OBJ_DIR)/IAR_calculate.o $(OBJ_DIR)/IAR_gradients.o $(OBJ_DIR)/Vector.o 
	
OBJS_ELLIPTIC=$(OBJ_DIR)/EBFV1_Elliptic_main.o $(OBJ_DIR)/EBFV1_AssemblyMatVec.o $(OBJ_DIR)/EBFV1_E.o $(OBJ_DIR)/EBFV1_G.o $(OBJ_DIR)/EBFV1_F_omega.o  $(OBJ_DIR)/EBFV1_wells.o \
              $(OBJ_DIR)/EBFV1_F_gamma.o $(OBJ_DIR)/EBFV1_MatrixFreeSolver.o $(OBJ_DIR)/EBFV1_PressureGradient.o $(OBJ_DIR)/MEBFV_Elliptic_main.o $(OBJ_DIR)/MEBFV_Initialize.o\
	          $(OBJ_DIR)/MEBFV_Assembly.o 

OBJS_HYPERBOLIC=$(OBJ_DIR)/EBFV1_velocityField.o $(OBJ_DIR)/EBFV1_hyperbolic.o $(OBJ_DIR)/EBFV1_hyperbolic_MIMPES.o $(OBJ_DIR)/EBFV1_advectiveTerm.o \
                $(OBJ_DIR)/EBFV1_advanceSaturation.o $(OBJ_DIR)/SaturationGradient.o $(OBJ_DIR)/EBFV1_hyperbolic_MIMPES_Adapt.o

OBJS_TOOLS=$(OBJ_DIR)/MeshData.o $(OBJ_DIR)/MeshData2.o $(OBJ_DIR)/OilProductionManagement.o $(OBJ_DIR)/LoadMeshAdaptationParameters.o  $(OBJ_DIR)/load_EBFV1_PreProcessorData.o $(OBJ_DIR)/auxiliar.o \
           $(OBJ_DIR)/SimulatorParameters.o $(OBJ_DIR)/SimulatorParameters2.o $(OBJ_DIR)/Exception.o $(OBJ_DIR)/exportVTK.o $(OBJ_DIR)/PhysicPropData.o $(OBJ_DIR)/PhysicPropData2.o $(OBJ_DIR)/LoadSimulatorParameters.o \
           $(OBJ_DIR)/Boundary_Conditions.o

OBJS_GEOMDATA=$(OBJ_DIR)/GeomData.o $(OBJ_DIR)/GeomData2.o $(OBJ_DIR)/GeomData3.o $(OBJ_DIR)/GeomData_initialize.o $(OBJ_DIR)/GeomData_calculate.o $(OBJ_DIR)/GeomData_Mapping.o \
              $(OBJ_DIR)/GeomData_AllocateMemory.o $(OBJ_DIR)/GeomData_TransferData.o 
	
OBJS_PREPROCESSOR=$(OBJ_DIR)/EBFV1__pre-processors.o $(OBJ_DIR)/EBFV1-2D-pp.o $(OBJ_DIR)/EBFV1-3D-pp.o $(OBJ_DIR)/EBFV1_modified.o $(OBJ_DIR)/rockProp.o \
	              $(OBJ_DIR)/assemblyMatrix_A.o $(OBJ_DIR)/calculateGeomCoefficients.o $(OBJ_DIR)/Matrix_E.o $(OBJ_DIR)/Matrix_G.o $(OBJ_DIR)/Matrix_F.o $(OBJ_DIR)/validate_coeff.o
	              
ifeq (,$(findstring NOADAPTATION,$(CXXFLAGS) $(ROCKPROP)))
OBJS_ADAPTATION=$(OBJ_DIR)/AdaptiveRemeshing.o \
    $(OBJ_DIR)/H_Refinement.o $(OBJ_DIR)/H_Refinement_2D.o \
    $(OBJ_DIR)/MeshRegularization.o $(OBJ_DIR)/SpecialElements_2D.o \
    $(OBJ_DIR)/UnevenElements.o $(OBJ_DIR)/FileManager.o $(OBJ_DIR)/FileManagerFunctions.o \
	$(OBJ_DIR)/MRE.o $(OBJ_DIR)/Rebuilder2.o $(OBJ_DIR)/ErrorEstimator.o $(OBJ_DIR)/Rebuilder3.o
endif

OBJECTS=$(OBJS_MAIN) $(OBJS_PREPROCESSOR) $(OBJS_ELLIPTIC) $(OBJS_HYPERBOLIC) $(OBJS_TOOLS) $(OBJS_GEOMDATA) $(OBJS_ADAPTATION) 

EXEC=PADMEC_AMR.exe
all:	$(EXEC)
#include ${PETSC_DIR}/bmake/common/base
include ${PETSC_DIR}/lib/petsc/conf/variables
include ${PETSC_DIR}/lib/petsc/conf/rules

$(EXEC):	$(OBJECTS) chkopts
	@echo "Linking objects..."
	-$(CXX) -o $(EXEC) $(OBJECTS) $(LIBS) $(PETSC_LIB)	

$(OBJ_DIR)/%.o:	$(SRC_DIR1)/%.cpp
	$(CXX) $(CXXFLAGS) $(ROCKPROP)  $(INCLUDES)  -c $<
	@mv *.o $(OBJ_DIR)
	
$(OBJ_DIR)/%.o:	$(SRC_DIR2)/%.cpp
	$(CXX) $(CXXFLAGS) $(ROCKPROP) $(INCLUDES)  -c $<
	@mv *.o $(OBJ_DIR)
	
$(OBJ_DIR)/%.o:	$(SRC_DIR21)/%.cpp
	$(CXX) $(CXXFLAGS) $(ROCKPROP) $(INCLUDES)  -c $< 
	@mv *.o $(OBJ_DIR)
	
$(OBJ_DIR)/%.o:	$(SRC_DIR4)/%.cpp
	$(CXX) $(CXXFLAGS) $(ROCKPROP) $(INCLUDES)  -c $< 
	@mv *.o $(OBJ_DIR)
	
$(OBJ_DIR)/%.o:	$(SRC_DIR5)/%.cpp
	$(CXX) $(CXXFLAGS) $(ROCKPROP) $(INCLUDES)  -c $< 
	@mv *.o $(OBJ_DIR)
	
$(OBJ_DIR)/%.o:	$(SRC_DIR51)/%.cpp
	$(CXX) $(CXXFLAGS) $(ROCKPROP) $(INCLUDES)  -c $< 
	@mv *.o $(OBJ_DIR)

$(OBJ_DIR)/%.o:	$(SRC_DIR6)/%.cpp
	$(CXX) $(CXXFLAGS) $(ROCKPROP) $(INCLUDES)  -c $< 
	@mv *.o $(OBJ_DIR)

$(OBJ_DIR)/%.o:	$(SRC_DIR7)/%.cpp
	$(CXX) $(CXXFLAGS) $(ROCKPROP) $(INCLUDES)  -c $< 
	@mv *.o $(OBJ_DIR)
	
$(OBJ_DIR)/%.o:	$(SRC_DIR71)/%.cpp
	$(CXX) $(CXXFLAGS) $(ROCKPROP) $(INCLUDES)  -c $< 
	@mv *.o $(OBJ_DIR)
	
$(OBJ_DIR)/%.o:	$(SRC_DIR8)/%.cpp
	$(CXX) $(CXXFLAGS) $(ROCKPROP) $(INCLUDES)  -c $< 
	@mv *.o $(OBJ_DIR)
	
$(OBJ_DIR)/%.o:	$(SRC_DIR10)/%.cpp
	$(CXX) $(CXXFLAGS) $(ROCKPROP) $(INCLUDES)  -c $< 
	@mv *.o $(OBJ_DIR)

ifeq (,$(findstring NOADAPTATION,$(CXXFLAGS) $(ROCKPROP)))
		
$(OBJ_DIR)/%.o:	$(SRC_DIR9)/%.cpp
	$(CXX) $(CXXFLAGS) $(ROCKPROP) $(INCLUDES)  -c $< 
	@mv *.o $(OBJ_DIR)

$(OBJ_DIR)/%.o:	$(SRC_DIR11)/%.cpp
	$(CXX) $(CXXFLAGS) $(ROCKPROP) $(INCLUDES)  -c $< 
	@mv *.o $(OBJ_DIR)

$(OBJ_DIR)/%.o:	$(SRC_DIR12)/%.cpp
	$(CXX) $(CXXFLAGS) $(ROCKPROP) $(INCLUDES)  -c $< 
	@mv *.o $(OBJ_DIR)
endif

#$(OBJ_DIR)/%.o:	$(SRC_DIR13)/%.cpp
#	$(CXX) $(CXXFLAGS) $(ROCKPROP) $(INCLUDES)  -c $< 
#	@mv *.o $(OBJ_DIR)
	
rebuild:
	@rm -f $(OBJ_DIR)/*.o $(EXEC)
	@echo "Limpeza concluida"
	@make
