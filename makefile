SRC:= $(wildcard src/*.cc)
OBJ:= $(SRC:src/%.cc=bin/%.o)

INCLUDE:= -Iinclude -IRLlib/include -Ilibconfig-1.4.9/lib
LIBS:= -lblas -lm -llapack -llapacke -LRLlib -lRLlib -Llibconfig-1.4.9/lib -Wl,-Bstatic -lconfig++ -Wl,-Bdynamic

MAINS := Compute
MAINSOBJ:= $(MAINS:%=bin/%.o)

GLOBALDEPEND:= $(wildcard include/*.hpp)

CCFLAGS:=  -pthread -Wall -fPIC -g -DHAVE_LAPACK_CONFIG_H -DLAPACK_COMPLEX_CPP

TEST:= test/RunTests

LIBCONFIG:= libconfig-1.4.9/lib/libconfig++.a
LIBCONFDIR:= libconfig-1.4.9

RLLIB:= RLlib/libRLlib.a

CC:= g++

.PHONY: clean doc all

all: bin $(RLLIB) $(LIBCONFIG) $(OBJ) $(MAINS) $(TEST)

$(LIBCONFIG): $(LIBCONFDIR)
	@echo Making libconfig...
	@(cd libconfig-1.4.9 && ./configure --prefix=$$(pwd) )
	@$(MAKE) -C libconfig-1.4.9
	@$(MAKE) -C libconfig-1.4.9 install

$(LIBCONFDIR): $(LIBCONFDIR).tar.gz
	@echo Unpacking libconfig...
	@tar -xzvf $<

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

$(RLLIB):
	$(MAKE) -C RLlib