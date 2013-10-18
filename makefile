SRC:= $(wildcard src/*.cc)
OBJ:= $(SRC:src/%.cc=bin/%.o)

INCLUDE:= -Iinclude -Ilapacke
LIBS:= -lblas -llapack -llapacke

MAINS := Compute
MAINSOBJ:= $(MAINS:%=bin/%.o)

GLOBALDEPEND:= $(wildcard include/*.hpp)

CCFLAGS:=  -pthread -Wall -fPIC -O3 -DHAVE_LAPACK_CONFIG_H -DLAPACK_COMPLEX_CPP

CC:= g++

.PHONY: clean doc


all: bin $(OBJ) $(MAINS)

bin: 
	mkdir -p bin

$(MAINS): %: $(filter-out $(MAINSOBJ), $(OBJ)) bin/%.o
	@echo Compiling $@...
	@$(CC) $(CCFLAGS) $(INCLUDE) $(LIBS) $^ -o $@  

$(OBJ): bin/%.o: src/%.cc include/%.hh $(GLOBALDEPEND)
	@echo Compiling $@...
	@$(CC) $(CCFLAGS) $(INCLUDE) $(LIBS) -c $< -o $@

clean:
	@echo Cleaning up...
	@rm -rf bin/*
	@rm -rf $(MAINSOBJ)
	@find . -type f -name "*~" -exec rm -f {} \;
	@find . -type f -name "\#*\#" -exec rm -f {} \;
