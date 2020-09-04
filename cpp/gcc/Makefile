RAD_CPP_DIR=	../
RAD_SRC_DIR=	$(RAD_CPP_DIR)/src
RAD_SRC_CORE_DIR=	$(RAD_SRC_DIR)/core
RAD_SRC_TRI_DIR=	$(RAD_SRC_DIR)/ext/triangle
SH_SRC_PARSE_DIR=	$(RAD_SRC_DIR)/ext/auxparse
SH_SRC_GEN_MATH_DIR=	$(RAD_SRC_DIR)/ext/genmath
RAD_SRC_ML_DIR=	$(RAD_SRC_DIR)/clients/mathlink
RAD_SRC_LIB_DIR=	$(RAD_SRC_DIR)/lib

#SHDIR=		$(SOFTDEVDIR)/Shared
#SHSRCDIR=	$(SHDIR)/src
#SHSRCGMDIR=	$(SHSRCDIR)/genmath
#SHSRCPARSEDIR=	$(SHSRCDIR)/auxparse
#SHLIBDIR=	$(SHDIR)/lib

LIB_DIR=	$(RAD_CPP_DIR)/../ext_lib

SH_INCLUDE=	-I$(SH_SRC_GEN_MATH_DIR) -I$(SH_SRC_PARSE_DIR)

MPI_INCLUDE_DIR_01=	/usr/lib/openmpi/include #NSLS-II cpu-001
MPI_INCLUDE_DIR_02=	/home/expdev/conda_envs/expdev_gcc7_openmp4/include #NSLS-II gpu-002
#Add more options here

MATHVERSION=	10.2
SYS=		Linux-x86-64

ML_ORIG_DIR=	/usr/local/Wolfram/Mathematica/${MATHVERSION}/SystemFiles/Links/MathLink/DeveloperKit/${SYS}/CompilerAdditions
#MMA_LIB_DIR=	/usr/local/Wolfram/Mathematica/${MATHVERSION}/SystemFiles/Libraries/${SYS}
MMA_LIB_DIR=	$(ML_ORIG_DIR)
#MLINCLUDEDIR=	$(ML_ORIG_DIR)
ML_INCLUDE=	-I$(ML_ORIG_DIR) -I$(RAD_SRC_ML_DIR)
#MLLIBDIR=	$(ML_ORIG_DIR)
#MLMPREP=	$(ML_ORIG_DIR)/mprep

#CC  = gcc
CC  = cc
#CC  = mpicxx
CXX = c++

#RAD_SRC_DEF=	-D__GCC__ -D__MATHEMATICA__ -DFFTW_ENABLE_FLOAT -DNO_TIMER -DANSI_DECLARATORS -DTRILIBRARY -DLINUX -D_GM_WITHOUT_BASE -DMLINTERFACE=4
#RAD_SRC_DEF=	-D__GCC__ -DFFTW_ENABLE_FLOAT -DNO_TIMER -DANSI_DECLARATORS -DTRILIBRARY -DLINUX -D_GM_WITHOUT_BASE -D_WITH_MPI
RAD_SRC_DEF=	-D__GCC__ -DFFTW_ENABLE_FLOAT -DNO_TIMER -DANSI_DECLARATORS -DTRILIBRARY -DLINUX -D_GM_WITHOUT_BASE

#CFLAGS=	-O3 -fPIC -I$(RAD_SRC_CORE_DIR) -I$(RAD_SRC_LIB_DIR) -I$(SH_SRC_PARSE_DIR) -I$(SH_SRC_GEN_MATH_DIR) -I$(RAD_SRC_TRI_DIR) $(ML_INCLUDE) $(RAD_SRC_DEF)
#LIBCFLAGS=	-O3 -fPIC -I$(RAD_SRC_CORE_DIR) -I$(RAD_SRC_LIB_DIR) -I$(SH_SRC_PARSE_DIR) -I$(SH_SRC_GEN_MATH_DIR) -I$(RAD_SRC_TRI_DIR) -I$(MPI_INCLUDE_DIR) $(RAD_SRC_DEF) -DALPHA__LIB__
LIBCFLAGS=	-O3 -fPIC -I$(RAD_SRC_CORE_DIR) -I$(RAD_SRC_LIB_DIR) -I$(SH_SRC_PARSE_DIR) -I$(SH_SRC_GEN_MATH_DIR) -I$(RAD_SRC_TRI_DIR) $(RAD_SRC_DEF) -DALPHA__LIB__

LIB_MATH=	-L$(MMA_LIB_DIR) -lML64i4 -ldl -luuid -L$(LIB_DIR) -lfftw -lm -lpthread -lrt
#LIB=		-L$(MMA_LIB_DIR) -l:/usr/local/Wolfram/Mathematica/10.2/SystemFiles/Links/MathLink/DeveloperKit/Linux-x86-64/CompilerAdditions/libML64i4.a -ldl -luuid -L$(LIB_DIR) -lfftw -lm -lpthread -lrt

LDFLAGS=	-L$(LIB_DIR) -lm -lfftw -lpthread -lrt -ldl -luuid

MODE ?= 0

ifeq ($(MODE), mpi)
LIBCFLAGS+= -I$(MPI_INCLUDE_DIR_01) -I$(MPI_INCLUDE_DIR_02) -D_WITH_MPI
else
ifeq ($(MODE), mpi_nersc)
LIBCFLAGS += -D_WITH_MPI
else
ifneq ($(MODE), 0)
$(error Unknown Radia compilation option)
endif
endif
endif

#OBJ=		radia_tm_linux.o radapl1.o radapl2.o radapl3.o radapl4.o radapl5.o radapl6.o radapl7.o radarccu.o radcast.o radexpgn.o radflm.o radg3d.o radg3dgr.o radgroup.o radinter.o radintrc.o radiobuf.o radmamet.o radmater.o radplnr1.o radplnr2.o radptrj.o radrec.o radrlmet.o radsbdac.o radsbdep.o radsbdrc.o radsbdvp.o radsend.o radvlpgn.o gmtrans.o gmfft.o triangle.o auxparse.o
OBJ=		radapl1.o radapl2.o radapl3.o radapl4.o radapl5.o radapl6.o radapl7.o radarccu.o radcast.o radexpgn.o radflm.o radg3d.o radg3dgr.o radgroup.o radinter.o radintrc.o radiobuf.o radmamet.o radmater.o radplnr1.o radplnr2.o radptrj.o radrec.o radrlmet.o radsbdac.o radsbdep.o radsbdrc.o radsbdvp.o radsend.o radvlpgn.o gmtrans.o gmfft.o triangle.o auxparse.o
OBJ_LIB=	radentry.o $(OBJ)
#OBJ_MATH=	radia_tm_linux.o $(OBJ)

PRG_LIB=		libradia.a
#PRG_MATH=		Radia.exe

lib:	$(OBJ_LIB)
	ar -cvq $(PRG_LIB) *.o
	#cp $(PRG_LIB) $(LIB_DIR)/
	rm -f *.o

%.o: $(RAD_SRC_LIB_DIR)/%.cpp
	$(CC) $(LIBCFLAGS) -c $<  

%.o: $(RAD_SRC_CORE_DIR)/%.cpp
	$(CC) $(LIBCFLAGS) -c $<  

triangle.o:	$(RAD_SRC_TRI_DIR)/triangle.c 
		cc $(LIBCFLAGS) -c $(RAD_SRC_TRI_DIR)/triangle.c -O $(RAD_SRC_DEF) -I$(RAD_SRC_TRI_DIR)

%.o: $(SH_SRC_PARSE_DIR)/%.cpp
	$(CC) $(LIBCFLAGS) -c $<  

%.o: $(SH_SRC_GEN_MATH_DIR)/%.cpp
	$(CC) $(LIBCFLAGS) -c $<  

#gmtrans.o:	$(SH_SRC_GEN_MATH_DIR)/gmtrans.cpp 
#		cc $(CFLAGS) -c $(SH_SRC_GEN_MATH_DIR)/gmtrans.cpp -O3 $(RAD_SRC_DEF)
#
#gmfft.o:	$(SH_SRC_GEN_MATH_DIR)/gmfft.cpp 
#		cc $(CFLAGS) -c $(SH_SRC_GEN_MATH_DIR)/gmfft.cpp -O3 $(RAD_SRC_DEF) 

clean:
	rm -f *.o *.so *.a

#pgm:	$(OBJ_MATH)
#	g++ $(ML_INCLUDE) $(OBJ_MATH) $(LIB_MATH) -o $(PRG_MATH)
#
#$(RAD_SRC_ML_DIR)/radia_tm_linux.c:	$(RAD_SRC_ML_DIR)/radia.tm
#		$(ML_ORIG_DIR)/mprep -o $(RAD_SRC_ML_DIR)/radia_tm_linux.c $(RAD_SRC_ML_DIR)/radia.tm
#
#radia_tm_linux.o:	$(RAD_SRC_ML_DIR)/radia_tm_linux.c 
#		cc $(CFLAGS) -c $(RAD_SRC_ML_DIR)/radia_tm_linux.c

clobber:
		@rm $(OBJ_MATH) $(PRG_MATH) radia_tm_linux.c *~ 

#all: lib pgm
all: lib

