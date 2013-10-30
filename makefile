SRC:= $(wildcard src/*.cc)
OBJ:= $(SRC:src/%.cc=bin/%.o)

INCLUDE:= -Iinclude -IRLlib/include
LIBS:= -lblas -lm -llapack -llapacke -LRLlib -lRLlib

MAINS := Compute
MAINSOBJ:= $(MAINS:%=bin/%.o)

GLOBALDEPEND:= $(wildcard include/*.hpp)

CCFLAGS:=  -pthread -Wall -fPIC -O3 -DHAVE_LAPACK_CONFIG_H -DLAPACK_COMPLEX_CPP

TEST:= test/RunTests

CC:= g++

.PHONY: clean doc

all: bin RLlib/libRLlib.a $(OBJ) $(MAINS) $(TEST)

bin: 
	mkdir -p bin

$(MAINS): %: $(filter-out $(MAINSOBJ), $(OBJ)) bin/%.o
	@echo Compiling $@...
	@$(CC) $(CCFLAGS) $(INCLUDE) $^ $(LIBS)  -o $@  

$(OBJ): bin/%.o: src/%.cc include/%.hh $(GLOBALDEPEND)
	@echo Compiling $@...
	@$(CC) $(CCFLAGS) $(INCLUDE) $(LIBS) -c $< -o $@

$(TEST): $(MAINS)
	@$(MAKE) -C test makerun

clean:
	@echo Cleaning up...
	@rm -rf bin/*
	@rm -rf $(MAINSOBJ)
	@find . -type f -name "*~" -exec rm -f {} \;
	@find . -type f -name "\#*\#" -exec rm -f {} \;
	@$(MAKE) -C RLlib clean
	@$(MAKE) -C test clean

RLlib/libRLlib.a:
	$(MAKE) -C RLlib