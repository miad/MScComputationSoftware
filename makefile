CC:= g++

CCFLAGS:=  -pthread -Wall -fPIC -DHAVE_LAPACK_CONFIG_H -DLAPACK_COMPLEX_CPP -DFP_SUPPORT_COMPLEX_DOUBLE_TYPE -DFP_USE_THREAD_SAFE_EVAL_WITH_ALLOCA -fstack-protector-all -std=c++11 -pedantic -g

CCOFLAGS:= -fexpensive-optimizations -fira-loop-pressure 

CCFLAGS+=$(CCOFLAGS)


GLOBALDEPEND:= $(wildcard include/*.hpp)


TEST:= test/RunTests

EXTRA_LIBS:= extra_libs

LIBCONFIG:= libconfig-1.4.9/lib/libconfig++.a
LIBCONFIGDIR:= libconfig-1.4.9

FPARSER:= fparser fpoptimizer fpconfig

FPARSERDIR:= fparser4.5.1

RLLIB:= RLlib/libRLlib.a


INCLUDE:= -Iinclude -IRLlib/include -I$(LIBCONFIGDIR)/lib
LIBS:= -lblas -lm -llapack -llapacke -LRLlib -lRLlib -Llibconfig-1.4.9/lib -Wl,-Bstatic -lconfig++ -Wl,-Bdynamic

SRC:= $(filter-out $(FPARSER:%=src/%.cc), $(wildcard src/*.cc)) $(FPARSER:%=src/%.cc)
OBJ:= $(SRC:src/%.cc=bin/%.o) 


MAINS := Compute
MAINSOBJ:= $(MAINS:%=bin/%.o)
MAINSNAME := Compute



.PHONY: clean doc all

all: bin $(RLLIB) $(LIBCONFIG) fparser_all $(OBJ) $(MAINS) $(TEST)


fparser_all: $(FPARSER:%=src/%.cc) $(FPARSER:%=include/%.hh) include/extrasrc

include/extrasrc: $(FPARSERDIR)/extrasrc
	@echo Linking $@...
	@ln -fs ../$< include/extrasrc

$(FPARSER:%=src/%.cc): $(FPARSERDIR)
	@echo Symlinking $@...
	@touch $(patsubst src/%, $</%, $@)
	@ln -fs $(patsubst src/%, ../$</%, $@) $@

$(FPARSER:%=include/%.hh): $(FPARSERDIR)
	@echo Symlinking $@...
	@touch $(patsubst include/%, $</%, $@)
	@ln -fs $(patsubst include/%, ../$</%, $@) $@


$(FPARSERDIR): $(EXTRA_LIBS)/$(FPARSERDIR).zip
	@echo Unpacking fparser...
	@mkdir -p $(FPARSERDIR)
	@unzip $< -d $@

$(LIBCONFIG): $(LIBCONFIGDIR)
	@echo Making libconfig...
	@(cd libconfig-1.4.9 && ./configure --prefix=$$(pwd) )
	@$(MAKE) -C libconfig-1.4.9
	@$(MAKE) -C libconfig-1.4.9 install

$(LIBCONFIGDIR): $(EXTRA_LIBS)/$(LIBCONFIGDIR).tar.gz
	@echo Unpacking libconfig...
	@tar -xzvf $<

bin: 
	mkdir -p bin

$(MAINS): %: $(filter-out $(MAINSOBJ), $(OBJ)) bin/%.o
	@echo Compiling $@...
	@$(CC) $(CCFLAGS) $(INCLUDE) $^ $(LIBS)  -o $(MAINSNAME)

$(OBJ): bin/%.o: src/%.cc include/%.hh $(GLOBALDEPEND)
	@echo Compiling $@...
	@$(CC) $(CCFLAGS) $(INCLUDE) $(LIBS) -c $< -o $@

$(TEST): $(MAINS)
	@$(MAKE) -C test makerun

clean: backup_clean fparse_clean
	@echo Cleaning up...
	@rm -rf bin/*
	@rm -rf $(MAINSOBJ)
	@$(MAKE) -C RLlib clean
	@$(MAKE) -C test clean

backup_clean: 
	@echo Removing tilde and hashtag files...
	@find . -type f -name "*~" -exec rm -f {} \;
	@find . -type f -name "\#*\#" -exec rm -f {} \;

fparse_clean:
	@echo Cleaning up fparser...
	@rm -f $(FPARSER:%=include/%.hh)
	@rm -f $(FPARSER:%=src/%.cc)
	@rm -rf include/extrasrc

$(RLLIB):
	$(MAKE) -C RLlib
