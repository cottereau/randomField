#### PRESUMED STRUCTURE
#		>>build        !(You are here to call the "make" command)
#		>>randomField  !(This is the git repository)
#			-BLAS
#			-SRC_EXE
#			-SRC_LIB
#		>>Tests !(Where the output will be written)
#

#### PATH TO HDF5 AND MPI LIBRARY AND INCLUDES - to be modified by user (to be completed)

#module load intel-compiler/15.0.1
#module load intel-mkl/11.2.1
#module load intel-mpi/5.0.2
#module load hdf5/1.8.12

#LIBHDF5 = -L/opt/san/bibliotheques/hdf5/1.8.12/lib/ -lhdf5 -lhdf5_hl -lhdf5_fortran
#INCLUDEHDF5 = -I/opt/san/bibliotheques/hdf5/1.8.12/include
#LIBMPI = -L/opt/san/intel/impi/4.0.0.028/lib64/ -lmpi -lmpi_dbg -lmpi_mt -lmpigf -lmpi_ilp64 
#INCLUDEMPI = -I/opt/san/intel/impi/4.0.0.028/include64 

LIBHDF5 = -L/opt/san/bibliotheques/hdf5/1.8.12/lib/ -lhdf5 -lhdf5_hl -lhdf5_fortran
INCLUDEHDF5 = -I/opt/san/bibliotheques/hdf5/1.8.12/include
LIBMPI = -L/opt/san/intel15/impi/5.0.2.044/lib64/ -lmpi -lmpi_dbg -lmpi_mt -lmpigf -lmpi_ilp64 
INCLUDEMPI = -I/opt/san/intel15/impi/5.0.2.044/include64

EXEC  = randomField.exe
#EXEC  = statistics.exe
#EXEC2 = statistics.exe
#FC   = ifort
FC   = mpiifort
FFLAGS = -O2

#### FROM THIS POINT ON THE USER SHOULD NOT MODIFY
RM := rm -rf

#Specifying the sources "f90" and "f"
SRCS = $(wildcard *.f90 *.f) 

#Create a ".f90" for each source
#OBJS = $(SRCS:.f90=.o) #SYNTAX NOT WORKING

OBJS += ./displayCarvalhol.o \
./main_RandomField.o \
./main_Statistics.o \
./test_func_RF.o \
./math_RF.o \
./constants_RF.o \
./randomFieldND.o \
./readFile_RF.o \
./spectra_RF.o \
./statistics_RF.o \
./writeResultFile_RF.o \
./lsame.o \
./charFunctions.o \
./xerbla.o \
./dger.o \
./mesh_RF.o \
./write_Log_File.o \
./systemUt_RF.o \
./common_variables_RF.o \
./type_RF.o \
./type_MESH.o \
./type_TEST.o \
./ranlib.o \
./rnglib.o \
./dgemm.o \
./ipmpar.o \
./spmpar.o \
./cumnor.o \
./devlpl.o \
./stvaln.o \
./dinvnr.o \
./cdfnor.o 

LIBS = $(LIBHDF5) $(LIBMPI) 
INCLUDE = $(INCLUDEHDF5) $(INCLUDEMPI)

#Dependencies
main_RandomField.o   : displayCarvalhol.o statistics_RF.o charFunctions.o \
					   randomFieldND.o  writeResultFile_RF.o readFile_RF.o \
					   test_func_RF.o constants_RF.o mesh_RF.o \
					   write_Log_File.o common_variables_RF.o \
					   systemUt_RF.o type_RF.o type_MESH.o type_TEST.o dger.o
main_Statistics.o    : displayCarvalhol.o statistics_RF.o charFunctions.o \
					   randomFieldND.o  writeResultFile_RF.o readFile_RF.o \
					   test_func_RF.o constants_RF.o mesh_RF.o \
					   write_Log_File.o common_variables_RF.o \
					   systemUt_RF.o type_RF.o type_MESH.o type_TEST.o dger.o					   
test_func_RF.o       : displayCarvalhol.o statistics_RF.o randomFieldND.o \
					   writeResultFile_RF.o readFile_RF.o mesh_RF.o mesh_RF.o\
					   write_Log_File.o common_variables_RF.o \
					   type_RF.o type_MESH.o type_TEST.o
writeResultFile_RF.o : displayCarvalhol.o math_RF.o statistics_RF.o write_Log_File.o constants_RF.o \
					   type_RF.o type_MESH.o type_TEST.o
randomFieldND.o      : displayCarvalhol.o math_RF.o spectra_RF.o \
			           constants_RF.o mesh_RF.o write_Log_File.o writeResultFile_RF.o\
			           common_variables_RF.o type_RF.o type_MESH.o dgemm.o ranlib.o
statistics_RF.o      : displayCarvalhol.o math_RF.o write_Log_File.o
spectra_RF.o         : displayCarvalhol.o math_RF.o write_Log_File.o constants_RF.o ranlib.o \
					   type_RF.o
mesh_RF.o            : math_RF.o write_Log_File.o type_RF.o type_MESH.o
math_RF.o            : displayCarvalhol.o write_Log_File.o constants_RF.o
systemUt_RF.o        : write_Log_File.o
write_Log_File.o     : charFunctions.o
type_RF.o            : charFunctions.o
type_MESH.o          : charFunctions.o
dgemm.o              : lsame.o xerbla.o
dger.o               : xerbla.o
ranlib.o             : rnglib.o

spmpar.o             : ipmpar.o
cumnor.o             : spmpar.o
stvaln.o             : devlpl.o
dinvnr.o             : stvaln.o cumnor.o
cdfnor.o             : dinvnr.o spmpar.o cumnor.o


# Making all the ".o" from the ".f90"
%.o: ../randomField/SRC_EXE/%.f90
	@echo 'Building file: $<'
	@echo 'Invoking: Fortran Compiler'
	$(FC) $(FFLAGS) -o "$@" $(INCLUDE) -c "$<"
	@echo 'Finished building: $<'
	@echo ' '
		
%.o: ../randomField/SRC_LIB/%.f90
	@echo 'Building file: $<'
	@echo 'Invoking: Fortran Compiler'
	$(FC) $(FFLAGS) -o "$@" $(INCLUDE) -c "$<"
	@echo 'Finished building: $<'
	@echo ' '

%.o: ../randomField/BLAS/%.f90
	@echo 'Building file: $<'
	@echo 'Invoking: Fortran Compiler'
	$(FC) $(FFLAGS) -o "$@" $(INCLUDE) -c "$<"
	@echo 'Finished building: $<'
	@echo ' '
	
%.o: ../randomField/RANDLIB/%.f90
	@echo 'Building file: $<'
	@echo 'Invoking: Fortran Compiler'
	$(FC) $(FFLAGS) -o "$@" $(INCLUDE) -c "$<"
	@echo 'Finished building: $<'
	@echo ' '	
	
%.o: ../randomField/dcdflib/src/%.f
	@echo 'Building file: $<'
	@echo 'Invoking: Fortran Compiler'
	$(FC) $(FFLAGS) -o "$@" $(INCLUDE) -c "$<"
	@echo 'Finished building: $<'
	@echo ' '	


	
# All Target
all: randomField

# Tool invocations
randomField: $(OBJS)
	@echo 'Building target: $@'
	@echo 'Invoking: Fortran Linker'
	$(FC) -o $(EXEC) $(FFLAGS) $(OBJS) $(INCLUDE) $(LIBS)
	@echo 'Finished building target: $@'
	@echo ' '	

# Other Targets
clean:
	-$(RM) $(EXECUTABLE) $(OBJS) *.mod
	-@echo ' '

.PHONY: all clean dependents
.SECONDARY:
