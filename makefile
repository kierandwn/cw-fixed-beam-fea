# MPI Compiler will compile all.
CC = mpicxx 
CXXFLAGS = -O3
BUILD = .build

INC = -I./includes
LDLIBS = -lscalapack-openmpi -lblacs-openmpi -lblacsCinit-openmpi -llapack -lblas -lboost_program_options
OBJS = ./$(BUILD)/$(EXEC).o ./$(BUILD)/assembly.o ./$(BUILD)/integration.o ./$(BUILD)/parallelise.o ./$(BUILD)/static.o

# All source / header files: not required by the file, included for user's benefit.
SRC  = assembly.cpp static.cpp integration.cpp parallelise.cpp
HDRS = assembly.h static.h integration.h parallelise.h
TXT  = static.txt dynamic.txt

# Name of the target executable.
EXEC = main

# Define standard arguments for program.
ARGS = --length 10000 --force 0.5 --area 12000.0 --I 14.4e+06 --E 210e+03 --T 3 --Tl 1 --rho 7850e-12

all: $(EXEC)

# help:

$(EXEC): $(OBJS) $(BUILD)
	$(CC) -o main $(OBJS) $(LOADLIBS) $(LDLIBS)

task1: $(EXEC)
	@ echo " > Computing Finite Element (Static) Solution... "
	./$(EXEC) $(ARGS) --static --Nx 24
	@ echo " > Passing to Python to compare with the analytical solution..."
	python ./post/static.py $(ARGS)
	@ echo " > Cleaning up intermediate files..."
	make -s clean

task2: $(EXEC)
	@ echo "Computing Finite Element (Dynamic) Solution: Explicit Time Integation Scheme... "
	./$(EXEC) $(ARGS) --dynamic --explicit --Nx 24 --Nt 20000
	@ echo "Processing the result with python..."
	python ./post/dynamic.py $(ARGS)
	@ echo "Cleaning up intermediate files..."
	make -s clean

task3: $(EXEC)
	@ echo " > Computing Finite Element (Dynamic) Solution: Implicit Time Integration Scheme... "
	./$(EXEC) $(ARGS) --dynamic --implicit --Nx 24 --Nt 20000
	@ echo " > Processing the result with python..."
	python ./post/dynamic.py $(ARGS)
	@ echo " > Cleaning up intermediate files..."
	make -s clean

task4: $(EXEC)
	@ echo " > Computing Finite Element (Dynamic) Solution: Explicit Time Integration Scheme (Parallelised)... "
	mpiexec -np 2 ./$(EXEC) $(ARGS) --dynamic --explicit --parallel --Nx 24 --Nt 100000
	@ echo " > Processing the result with python..."
	python ./post/dynamic.py $(ARGS)
	@ echo " > Cleaning up intermediate files..."
	make -s clean

task5: $(EXEC)
	@ echo " > Computing Finite Element (Dynamic) Solution: Implicit Time Integration Scheme (Parallelised)... "
	mpiexec -np 2 ./$(EXEC) $(ARGS) --dynamic --implicit --parallel --Nx 24 --Nt 100000
	@ echo " > Processing the result with python..."
	python ./post/dynamic.py $(ARGS)
	@ echo " > Cleaning up intermediate files..."
	make -s clean

# Target specific to main files since it is stored in root directory.
./$(BUILD)/$(EXEC).o : $(EXEC).cpp dir
	@ echo " > Compiling... "
	$(CC) $(CXXFLAGS) $(INC) -o $@ -c $<

./$(BUILD)/%.o : ./src/%.cpp dir
	$(CC) $(CXXFLAGS) $(INC) -o $@ -c $<

.PHONY: dir
dir: $(BUILD)

$(BUILD):
	@ echo " > Creating build folder..."
	mkdir -p $(BUILD)

.PHONY: clean
clean:
	rm -r $(BUILD)

.PHONY: cleanall
cleanall:
	-rm -r -f $(BUILD) $(EXEC)
	rm -f ./output/*.txt