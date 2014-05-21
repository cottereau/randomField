################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
F90_SRCS += \
../displayCarvalhol.f90 \
../mainCA1D.f90 \
../randField1D.f90 \
../readRF.f90 \
../vectorOpCarvalhol.f90 

OBJS += \
./displayCarvalhol.o \
./mainCA1D.o \
./randField1D.o \
./readRF.o \
./vectorOpCarvalhol.o 


# Each subdirectory must supply rules for building sources it contributes
%.o: ../%.f90
	@echo 'Building file: $<'
	@echo 'Invoking: GNU Fortran Compiler'
	gfortran -funderscoring -I/opt/san/bibliotheques/hdf5/1.8.12/include -I/opt/san/compilateurs/gcc-4.8.2/include -I/opt/san/bibliotheques/mpc/1.0.2/gcc/include -I/opt/san/bibliotheques/mpfr/3.1.2/gcc/include -I/opt/san/intel13/2013.0.028/composer_xe_2013_sp1/mkl/include -I/opt/san/intel/impi/4.0.0.028/include64 -O0 -g -Wall -c -fmessage-length=0 -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '

displayCarvalhol.o: ../displayCarvalhol.f90

mainCA1D.o: ../mainCA1D.f90 displayCarvalhol.o randField1D.o readRF.o vectorOpCarvalhol.o

randField1D.o: ../randField1D.f90 displayCarvalhol.o vectorOpCarvalhol.o

readRF.o: ../readRF.f90 displayCarvalhol.o

vectorOpCarvalhol.o: ../vectorOpCarvalhol.f90


