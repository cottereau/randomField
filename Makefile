################################################################################
# Automatically-generated file. Do not edit!
################################################################################

-include ../makefile.init

RM := rm -rf

# All of the sources participating in the build are defined here

F90_SRCS += \
../displayCarvalhol.f90 \
../main_RandomField.f90 \
../randomField1D.f90 \
../randomField2D.f90 \
../randomField3D.f90 \
../readInputFile_RF.f90 \
../spectra_RF.f90 \
../statistics_RF.f90 \
../writeResultFile_RF.f90 

OBJS += \
./displayCarvalhol.o \
./main_RandomField.o \
./randomField1D.o \
./randomField2D.o \
./randomField3D.o \
./readInputFile_RF.o \
./spectra_RF.o \
./statistics_RF.o \
./writeResultFile_RF.o 


# Each subdirectory must supply rules for building sources it contributes
%.o: ../%.f90
	@echo 'Building file: $<'
	@echo 'Invoking: GNU Fortran Compiler'
	gfortran -funderscoring -O0 -g -Wall -c -fmessage-length=0 -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '

displayCarvalhol.o: ../displayCarvalhol.f90

main_RandomField.o: ../main_RandomField.f90 randomField1D.o statistics_RF.o

randomField1D.o: ../randomField1D.f90 displayCarvalhol.o spectra_RF.o

randomField2D.o: ../randomField2D.f90

randomField3D.o: ../randomField3D.f90

readInputFile_RF.o: ../readInputFile_RF.f90

spectra_RF.o: ../spectra_RF.f90

statistics_RF.o: ../statistics_RF.f90

writeResultFile_RF.o: ../writeResultFile_RF.f90


-include ../makefile.defs

# Add inputs and outputs from these tool invocations to the build variables 

# All Target
all: randomField

# Tool invocations
randomField: $(OBJS) $(USER_OBJS)
	@echo 'Building target: $@'
	@echo 'Invoking: GNU Fortran Linker'
	gfortran  -o "randomField" $(OBJS) $(USER_OBJS) $(LIBS)
	@echo 'Finished building target: $@'
	@echo ' '

# Other Targets
clean:
	-$(RM) $(EXECUTABLES)$(OBJS)$(C_DEPS) randomField
	-@echo ' '

.PHONY: all clean dependents
.SECONDARY:

-include ../makefile.targets
