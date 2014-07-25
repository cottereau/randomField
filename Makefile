#### PRESUMED STRUCTURE
#		>>Build
#		>>randomField  !(This is the git repository)
#			-Source
#		>>Tests !(Where the output will be written)
#

#### PATH TO HDF5 AND MPI LIBRARY AND INCLUDES - to be modified by user (to be completed)

#module load intel-compiler/14.0.0
#module load intel-mkl/11.1.0
#module load intel-mpi/4.0.0.028
#module load hdf5/1.8.12

LIBHDF5 = -L/opt/san/bibliotheques/hdf5/1.8.12/lib/ -lhdf5 -lhdf5_hl -lhdf5_fortran
INCLUDEHDF5 = -I/opt/san/bibliotheques/hdf5/1.8.12/include
LIBMPI = -L/opt/san/intel/impi/4.0.0.028/lib64/ -lmpi -lmpi_dbg -lmpi_mt -lmpigf -lmpi_ilp64 
INCLUDEMPI = -I/opt/san/intel/impi/4.0.0.028/include64 

EXEC = randomField.exe
#FC   = ifort
FC   = mpiifort

#### FROM THIS POINT ON THE USER SHOULD NOT MODIFY
RM := rm -rf

#Specifying the sources "f90"
SRCS = $(wildcard *.f90)
#./displayCarvalhol.f90 \
#./main_RandomField.f90 \
#./math_RF.f90 \
#./randomFieldND.f90 \
#./readInputFile_RF.f90 \
#./spectra_RF.f90 \
#./statistics_RF.f90 \
#./writeResultFile_RF.f90

#Create a ".f90" for each source
#OBJS = $(SRCS:.f90=.o) #SYNTAX NOT WORKING
OBJS += ./displayCarvalhol.o \
./main_RandomField.o \
./math_RF.o \
./randomFieldND.o \
./readFile_RF.o \
./spectra_RF.o \
./statistics_RF.o \
./writeResultFile_RF.o

LIBS = $(LIBHDF5) $(LIBMPI) 
INCLUDE = $(INCLUDEHDF5) $(INCLUDEMPI)

# Making all the ".o" from the ".f90"
%.o: ../randomField/%.f90
	@echo 'Building file: $<'
	@echo 'Invoking: Fortran Compiler'
	$(FC) -o "$@" $(INCLUDE) -c "$<"
	@echo 'Finished building: $<'
	@echo ' '

#Dependencies
main_RandomField.o   : statistics_RF.o randomFieldND.o  writeResultFile_RF.o readFile_RF.o
writeResultFile_RF.o : displayCarvalhol.o math_RF.o statistics_RF.o
randomFieldND.o      : displayCarvalhol.o math_RF.o spectra_RF.o    
statistics_RF.o      : displayCarvalhol.o math_RF.o
spectra_RF.o         : displayCarvalhol.o math_RF.o

# All Target
all: randomField

# Tool invocations
randomField: $(OBJS)
	@echo 'Building target: $@'
	@echo 'Invoking: Fortran Linker'
	$(FC) -o $(EXEC) $(OBJS) $(INCLUDE) $(LIBS)
	@echo 'Finished building target: $@'
	@echo ' '

# Other Targets
clean:
	-$(RM) $(EXECUTABLE) $(OBJS) *.mod
	-@echo ' '

.PHONY: all clean dependents
.SECONDARY:
