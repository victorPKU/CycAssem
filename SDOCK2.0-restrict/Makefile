#                              SDOCK2.0
#   Global Protein-Protein Docking Using Stepwise Force-Field Potentials   
#                Written by ZHANG Changsheng 
#  Copyright (c) 2024 Molecular Design Lab at Peking University, Beijing, China  
 
#    contact: changshengzhang@pku.edu.cn, lhlai@pku.edu.cn

# please change the FFTW library route
FFTW_DIR        = /home/changsheng/Downloads/Softwares/MDL/SDOCK1p0/fftw-3.3.8


SHELL           = /bin/sh

CC              = gcc

CC_FLAGS        = -ansi -O6 -fomit-frame-pointer -Wall -W -Wcast-qual -Wpointer-arith -Wcast-align -pedantic -fno-schedule-insns -fschedule-insns2 -fstrict-aliasing

CC_LINKERS      = -lm

STRIP           = strip

SECURITY	= chmod 555


CC_FLAGS_FULL	= -I$(FFTW_DIR)/include $(CC_FLAGS)
FFTW_LINKERS    = -L$(FFTW_DIR)/.libs -lfftw3f 


.SUFFIXES:	.c .o

.c.o:
		$(CC) $(CC_FLAGS_FULL) -c $<


PROGRAMS = preprocess watmap sdock build

all:		$(PROGRAMS)
 

preprocess:		preprocess.o protein.o structureIO.o surface.o mem.o GNM.o geometry.o mem.h protein.h surface.h GNM.h geometry.h
		$(CC) $(CC_FLAGS) -o $@ preprocess.o protein.o structureIO.o surface.o mem.o GNM.o geometry.o $(CC_LINKERS)
		$(STRIP) $@
		$(SECURITY) $@


watmap:		watmap.o watgrid.o structure.o mem.o geometry.o watgrid.h mem.h structure.h geometry.h
		$(CC) $(CC_FLAGS) -o $@ watmap.o watgrid.o structure.o mem.o geometry.o $(CC_LINKERS)
		$(STRIP) $@
		$(SECURITY) $@


sdock:		sdock.o cluster.o fftmem.o grid.o mem.o record.o rotplan.o structure.o geometry.o cluster.h fftmem.h grid.h mem.h record.h rotplan.h structure.h geometry.h
		$(CC) $(CC_FLAGS_FULL) -o $@ sdock.o cluster.o fftmem.o grid.o mem.o record.o rotplan.o structure.o geometry.o $(FFTW_LINKERS) $(CC_LINKERS)
		$(STRIP) $@
		$(SECURITY) $@


build:		build.o mem.o model.o record.o rotplan.o structure.o geometry.o build.h mem.h record.h rotplan.h structure.h geometry.h
		$(CC) $(CC_FLAGS) -o $@ build.o mem.o model.o record.o rotplan.o structure.o geometry.o $(LIBRARY_OBJECTS) $(CC_LINKERS)
		$(STRIP) $@
		$(SECURITY) $@

clean:
		rm -f *.o core $(PROGRAMS)


preprocess.o:		protein.h mem.h
protein.o:			protein.h mem.h surface.h
structureIO.o:		protein.h mem.h
surface.o:			surface.h mem.h
GNM.o:			GNM.h mem.h geometry.h protein.h
mem.o:			mem.h
geometry.o:              geometry.h
watmap.o:			watgrid.h structure.h geometry.h mem.h
structure.o:		structure.h mem.h 
watgrid.o:			watgrid.h structure.h geometry.h mem.h
sdock.o:			grid.h structure.h record.h cluster.h rotplan.h fftmem.h mem.h
grid.o:			grid.h structure.h
cluster.o:			cluster.h structure.h rotplan.h record.h mem.h
record.o:			record.h mem.h
fftmem.o:			fftmem.h
rotplan.o:			rotplan.h mem.h
build.o:			build.h structure.h rotplan.h record.h mem.h
model.o:            build.h structure.h

