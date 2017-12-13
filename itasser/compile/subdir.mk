################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
O_SRCS += \
../cas_mod.o 

F90_SRCS += \
../cas_mod.f90 \
../nmr.f90 

OBJS += \
./cas_mod.o \
./nmr.o 


# Each subdirectory must supply rules for building sources it contributes
%.o: ../%.f90
	@echo 'Building file: $<'
	@echo 'Invoking: GNU Fortran Compiler'
	gfortran -funderscoring -O3 -Wall -c -fmessage-length=0 -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '

cas_mod.o: ../cas_mod.f90 nmr.o

nmr.o: ../nmr.f90


