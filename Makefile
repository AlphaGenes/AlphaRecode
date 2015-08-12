# General variables
NAME:=AlphaRecode
VERSION:=1
SUBVERSION:=0
PROGRAM:=$(NAME)$(VERSION).$(SUBVERSION)

# Set the default compiler to iFort
FC:=ifort
FFLAGS:=-O3

# Flags required for C++ Formatter compilation
CFLAGS:= -std=c++0x -Wextra -Wall 
CC:=g++

ifeq ($(OS), Windows_NT)
	OSFLAG := "OS_WIN"
	FFLAGS := $(FFLAGS) /static /i8 /fpp /Qopenmp-link:static /module $(BUILDDIR) /Qmkl /Qlocation,link,"$(VCINSTALLDIR)/bin" /D $(OSFLAG)
	ABOPT := -static -Qopenmp-link:static -Qmkl -Qlocation,link,"$(VCINSTALLDIR)/bin"
	obj := .obj
	exe := .exe
else
	# Linux or Mac OSX
	obj := .o
	OSFLAG := "OS_UNIX"
	ABOPT := -mkl -static-intel -openmp-link=static
	exe :=
	FFLAGS:= $(FFLAGS) -mkl -i8 -static-intel -fpp -openmp-link=static  -D $(OSFLAG)
	uname := $(shell uname)
  # Linux only
	ifeq ($(uname), Linux)
		FFLAGS := $(FFLAGS) -static -static-libgcc -static-libstdc++
	endif
endif

all: $(NAME)$(exe)


debug: FFLAGS:= -i8 -traceback -g -D $(OSFLAG) -debug all -warn -check bounds -check format \
		-check output_conversion -check pointers -check uninit -fpp

debug: all 

$(NAME)$(exe): $(SRCDIR)$(NAME).f90
	@echo "Compiling $(NAME)..."
	$(FC)  $(SRCDIR)$(NAME).f90 $(FFLAGS) -o $(NAME)$(exe)
	@echo

clean:
	rm -rf *$(obj) *.mod *.dwarf *.i90 *__genmod* *~ $(NAME)$(exe)