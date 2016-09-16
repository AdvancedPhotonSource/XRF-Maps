USEBLAS = -D_HAVE_CBLAS
#LIBS = -L/kyb/agbs/suvrit/work/src/GotoBLAS2 -lgoto2 -lgfortran
# If you are using a blas library, then include the following, else use
# LIBS = 
LIBS = -lblas
MACHINE = -m64
#MACHINE = -m32

CXXFLAGS = -Wall -pg -Wno-write-strings -funroll-loops -O2 $(MACHINE) -D_LARGEFILE64_SOURCE  -D_FILE_OFFSET_BITS=64 $(USEBLAS)
CXX = g++

# If your cblas.h file resides in a strange place
INC = 

%.o: %.cc %.h
	$(CXX) -c $(INC) $(CXXFLAGS)  $< -o $@

LIBOBJS = denseMatrix.o sparseMatrix.o nnls.o

all: libnnls.a nnls

libnnls.a: $(LIBOBJS)
	ar cr libnnls.a $(LIBOBJS) 
	ranlib libnnls.a

nnls: libnnls.a nnlsDriver.o
	$(CXX) $(INC) $(CXXFLAGS) nnlsDriver.o -o nnls -L. -lnnls $(LIBS)

clean:
	rm -f *.o libnnls.a


VER = 0.52
FILES = nnls.cc nnls.h nnlsDriver.cc denseMatrix.h sparseMatrix.h matrix.h vector.h makefile denseMatrix.cc sparseMatrix.cc Copyright LICENSE README pbbnnls.m solopt.m
UTIL = loadCmatrix.m matlab2ccs.c readCCS.c saveMatlabMatrix.m
EGS = tst_col_ccs tst_dim tst_row_ccs tst.solution tst_txx_nz
DEST = nnls-$(VER)
mvhere = $(PWD)

dist:
	rm -rf /tmp/$(DEST)
	mkdir -p /tmp/$(DEST)
	mkdir -p /tmp/$(DEST)/util
	mkdir -p /tmp/$(DEST)/examples
	cp $(FILES) /tmp/$(DEST)
	cd util && cp $(UTIL) /tmp/$(DEST)/util/
	cd examples && cp $(EGS) /tmp/$(DEST)/examples
	cd /tmp && tar cvzf $(DEST).tar.gz $(DEST) && mv $(DEST).tar.gz $(mvhere)
