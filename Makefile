##  file:    Makefile
##  update:  09/24/03

default: esolid

ESOLID_dir = .

###  1. System configuration

####  1-1. Hardware and operating system

#####  Pick your architecture.
CFLAGS = -D_Linux_i386_
#CFLAGS = -D_EXPERIMENT -D_Linux_i386_
#CFLAGS = -D_SunOS_
#CFLAGS = -D_IRIX_

####  1-2. Compiler

#####  Pick your compiler.
#####  Currently, the library is tested only with g++.
CC = g++
#CC = g++-2.95
#CC = g++-3.0
#CC = g++-3.2

#####  Define g++ compiler flags.
CFLAGS += -g
CFLAGS += -O
#CFLAGS += -O2
#CFLAGS += -Wno-deprecated

#####  Pick your archiver and stripper.

######  "ar -ru" replaces the objects only if they are updated.
#AR = ar -r
AR = ar -ru

######  "strip --strip-unneeded" strips all symbols
######    that are not needed for relocation.
######  "strip -g" strips debugging symbols only.
#STRIP = strip --strip-unneeded
STRIP = strip -g --strip-unneeded
######  For debugging, do not strip, so use the "echo" (dummy) definition: 
#STRIP = echo


###  3. Header files and library locations

####  3-1. GMP

#IFLAGS_gmp = -I/usr/include
LFLAGS_gmp = -lgmp

####  3-2. LPACK

#IFLAGS_lapack = -I/usr/local/include
#LFLAGS_lapack = -L/usr/local/lib -lclapack -lcblas -lF77
LFLAGS_lapack = -lclapack -lcblas -lF77

###  4. Making

####  4-1. Objects

%.o: %.cc
	$(CC) $(CFLAGS) -c $(ESOLID_dir)/$< -o $(ESOLID_dir)/$(basename $<).o

config_objects = \
	$(ESOLID_dir)/src/timer.o \
	$(ESOLID_dir)/src/config.o \

bignum_objects = \
	$(ESOLID_dir)/src/bigint.o \
	$(ESOLID_dir)/src/bigint_vector.o \
	$(ESOLID_dir)/src/bigint_matrix.o \
	$(ESOLID_dir)/src/bigrational.o \
	$(ESOLID_dir)/src/bigrational_vector.o \
	$(ESOLID_dir)/src/bigrational_matrix.o \

fp_objects = \
	$(ESOLID_dir)/src/fpfilter.o \
	$(ESOLID_dir)/src/fpconversion.o \

mapc_objects = \
	$(ESOLID_dir)/src/kratpoly.o \
	$(ESOLID_dir)/src/kfloatpoly.o \
	$(ESOLID_dir)/src/root1.o \
	$(ESOLID_dir)/src/VDM.o \
	$(ESOLID_dir)/src/kpoint1d.o \
	$(ESOLID_dir)/src/kpoint2d.o \
	$(ESOLID_dir)/src/kboxco2.o \
	$(ESOLID_dir)/src/ksegment.o \
	$(ESOLID_dir)/src/kcurve.o \

esolid_objects = \
	$(ESOLID_dir)/src/kbox3d.o \
	$(ESOLID_dir)/src/pt_surf_assoc.o \
	$(ESOLID_dir)/src/pascal.o \
	$(ESOLID_dir)/src/ksurf.o \
	$(ESOLID_dir)/src/kpatch.o \
	$(ESOLID_dir)/src/kpartition.o \
	$(ESOLID_dir)/src/kgraph.o \
	$(ESOLID_dir)/src/genbox.o \
	$(ESOLID_dir)/src/gencyl.o \
	$(ESOLID_dir)/src/genell.o \
	$(ESOLID_dir)/src/gentor.o \
	$(ESOLID_dir)/src/ksolid.o \

####  4-2. Archives

IFLAGS  = -Iinclude $(IFLAGS_gmp) $(IFLAGS_lapack)
CFLAGS += $(IFLAGS)

LFLAGS  = $(LFLAGS_gmp) $(LFLAGS_lapack)
LFLAGS += -lm

#####  Make mapc static library.
libmapc:	$(config_objects) $(bignum_objects) $(fp_objects) $(mapc_objects)
	$(AR) $(ESOLID_dir)/$@.a $(config_objects) $(bignum_objects) $(fp_objects) $(mapc_objects)
	$(STRIP) $(ESOLID_dir)/$@.a

#####  Make mapc shared library.
#libmapc:	$(config_objects) $(bignum_objects) $(fp_objects) $(mapc_objects)
#	$(CC) -shared -Wl -o $(ESOLID_dir)/$@.so $(config_objects) $(bignum_objects) $(fp_objects) $(mapc_objects)
#	$(STRIP) $(ESOLID_dir)/$@.so

#####  Make esolid static library.
libesolid:	$(config_objects) $(bignum_objects) $(fp_objects) $(mapc_objects) $(esolid_objects)
	$(AR) $(ESOLID_dir)/$@.a $(config_objects) $(bignum_objects) $(fp_objects) $(mapc_objects) $(esolid_objects)
	$(STRIP) $(ESOLID_dir)/$@.a

#####  Make esolid shared library
#libesolid:	$(config_objects) $(bignum_objects) $(fp_objects) $(mapc_objects) $(esolid_objects)
#	$(CC) -shared -Wl -o $(ESOLID_dir)/$@.so $(config_objects) $(bignum_objects) $(fp_objects) $(mapc_objects) $(esolid_objects)
#	$(STRIP) $(ESOLID_dir)/$@.so

####  4-3.  main functions

#####  Make executable

mapc_poly:	$(ESOLID_dir)/main/mapc_poly_main.o libmapc
	$(CC) -o $(ESOLID_dir)/$@ $(ESOLID_dir)/main/mapc_poly_main.o -L$(ESOLID_dir) -lmapc  $(CFLAGS) $(LFLAGS)

mapc_sturm:	$(ESOLID_dir)/main/mapc_sturm_main.o libmapc
	$(CC) -o $(ESOLID_dir)/$@ $(ESOLID_dir)/main/mapc_sturm_main.o -L$(ESOLID_dir) -lmapc $(CFLAGS) $(LFLAGS)

#mapc_rt1:	$(ESOLID_dir)/main/mapc_rt1_main.o libmapc
#	$(CC) -o $(ESOLID_dir)/$@ $(ESOLID_dir)/main/mapc_rt1_main.o -L$(ESOLID_dir) -lmapc  $(CFLAGS) $(LFLAGS)

mapc_pt1:	$(ESOLID_dir)/main/mapc_pt1_main.o libmapc
	$(CC) -o $(ESOLID_dir)/$@ $(ESOLID_dir)/main/mapc_pt1_main.o -L$(ESOLID_dir) -lmapc  $(CFLAGS) $(LFLAGS)

#mapc_pt2:	$(ESOLID_dir)/main/mapc_pt2_main.o libmapc
#	$(CC) -o $(ESOLID_dir)/$@ $(ESOLID_dir)/main/mapc_pt2_main.o -L$(ESOLID_dir) -lmapc  $(CFLAGS) $(LFLAGS)

mapc:		$(ESOLID_dir)/main/mapc_main.o libmapc
	$(CC) -o $(ESOLID_dir)/$@ $(ESOLID_dir)/main/mapc_main.o -L$(ESOLID_dir) -lmapc  $(CFLAGS) $(LFLAGS)

esolid:		$(ESOLID_dir)/main/esolid_main.o libesolid
	$(CC) -o $(ESOLID_dir)/$@ $(ESOLID_dir)/main/esolid_main.o -L$(ESOLID_dir) -l$@  $(CFLAGS) $(LFLAGS)

###  5. Cleaning

clean:
	rm -f $(ESOLID_dir)/src/*.o $(bignum_objects) $(esolid_objects) $(general_objects) $(mapc_objects) $(pseudoroot_objects)
	rm -f $(ESOLID_dir)/libmapc.a $(ESOLID_dir)/libmapc.so $(ESOLID_dir)/libesolid.a $(ESOLID_dir)/libesolid.so
	rm -f $(ESOLID_dir)/main/*.o
	rm -f $(ESOLID_dir)/mapc* $(ESOLID_dir)/esolid

