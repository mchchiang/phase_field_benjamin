
CC := gcc # GNU compiler
# CC := icc # Intel compiler 
CFLAGS := -Wall -Wshadow -Winline -O3 -std=c11 # No OpenMP
#CFLAGS := -Wall -Wshadow -Winline -O3 -std=c11 -fopenmp # GNU OpenMP
#CFLAGS := -Wall -Wshadow -Winline -O3 -std=c11 -qopenmp # Intel OpenMP


CXX := g++ # GNU compiler
# CXX := icpc # Intel compiler
CXXFLAGS := -Wall -Wshadow -Winline -O3 -std=c++11 # No OpenMP
#CXXFLAGS := -Wall -Wshadow -Winline -O3 -std=c++11 -fopenmp # GNU OpenMP
#CXXFLAGS := -Wall -Wshadow -Winline -O3 -std=c++11 -qopenmp # Intel OpenMP

# 0 = no armadillo and 1 = with armadillo
HAS_ARMA := 0

INC := 
LFLAGS := -lm
#LFLAGS := -lm -lgomp

PRG_SUFFIX_FLAG := 0
PRG_SUFFIX =

BDIR := ../bin
ODIR := $(BDIR)/obj/
EDIR := $(BDIR)/exe/

SRC_C := arralloc.c mtwister.c cell.c phase_field_model.c array.c image.c overlap.c dump.c dump_cm.c dump_bulk_cm.c dump_field.c dump_cell_field.c dump_gyration.c dump_neighbour.c dump_energy.c dump_overlap.c dump_overlap_field.c dump_index_field.c random_mtwister.c
SRC_CXX := neighbour.cpp
ifeq ($(HAS_ARMA),1)
	SRC_C := $(SRC_C) shape.c dump_shape.c
	SRC_CXX := $(SRC_CXX) sgolay_filter.cpp
	LFLAGS := $(LFLAGS) -larmadillo
endif
OBJ_C := $(patsubst %.c, $(ODIR)/%.o, $(SRC_C))
OBJ_CXX := $(patsubst %.cpp, $(ODIR)/%.o, $(SRC_CXX))
PRG_SRC := run_phase_field_model.c
PRG_OBJ := $(patsubst %.c, $(ODIR)/%.o, $(PRG_SRC))
PRG := $(patsubst %.c, $(EDIR)/%$(PRG_SUFFIX), $(PRG_SRC))

all : check_arma make_bin $(LOBJ) $(OBJ) $(PRG_OBJ) $(PRG)

%.o : %.c %.h
	$(CC) $(CFLAGS) $(INC) -o $@ -c $<

$(ODIR)/%.o : %.c %.h
	$(CC) $(CFLAGS) $(INC) -o $@ -c $<

$(ODIR)/%.o : %.c
	$(CC) $(CFLAGS) $(INC) -o $@ -c $<

$(ODIR)/%.o : %.cpp %.hpp
	$(CXX) $(CXXFLAGS) $(INC) -o $@ -c $<

$(ODIR)/%.o : %.cpp %.h
	$(CXX) $(CXXFLAGS) $(INC) -o $@ -c $<

$(ODIR)/%.o : %.cpp
	$(CXX) $(CXXFLAGS) $(INC) -o $@ -c $<

$(PRG) : $(OBJ_C) $(OBJ_CXX) $(PRG_OBJ)
	$(CXX) $(LFLAGS)  $^ -o $@

.PHONY : clean check_arma make_bin

clean :
	rm -rf $(OBJ_C) $(OBJ_CXX) $(PRG_OBJ) $(PRG)
	rm -rf $(ODIR) $(EDIR) $(BDIR)

check_arma :
ifeq ($(HAS_ARMA),1)
	perl -pi -e 's/PF_HAS_ARMA [0,1]/PF_HAS_ARMA 1/g' dump.h
else
	perl -pi -e 's/PF_HAS_ARMA [0,1]/PF_HAS_ARMA 0/g' dump.h
endif

make_bin :
	mkdir -p $(ODIR)
	mkdir -p $(EDIR)
