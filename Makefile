#### PATH TO HDF5 LIBRARY - to be modified by user (to be completed)
LIBHDF5 = 

EXECUTABLE = randomField.x

RM := rm -rf

OBJS += displayCarvalhol.o \
./main_RandomField.o \
./math_RF.o \
./randomFieldND.o \
./readInputFile_RF.o \
./spectra_RF.o \
./statistics_RF.o \
./writeResultFile_RF.o

LIBS = $(LIBHDF5) 

# Each subdirectory must supply rules for building sources it contributes
%.o: ../randomField/%.f90
	@echo 'Building file: $<'
	@echo 'Invoking: GNU Fortran Compiler'
	gfortran -funderscoring -O0 -g -Wall -c -fmessage-length=0 -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '

main_RandomField.o: randomFieldND.o statistics_RF.o writeResultFile_RF.o

randomFieldND.o: displayCarvalhol.o spectra_RF.o math_RF.o

spectra_RF.o: displayCarvalhol.o math_RF.o

# All Target
all: randomField

# Tool invocations
randomField: $(OBJS)
	@echo 'Building target: $@'
	@echo 'Invoking: GNU Fortran Linker'
	gfortran -o $(EXECUTABLE) $(OBJS) $(LIBS)
	@echo 'Finished building target: $@'
	@echo ' '

# Other Targets
clean:
	-$(RM) $(EXECUTABLE) $(OBJS) *.mod
	-@echo ' '

.PHONY: all clean dependents
.SECONDARY:
