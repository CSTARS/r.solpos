MODULE_TOPDIR = ../..

PGM = r.solpos

LIBES = $(GPROJLIB) $(RASTERLIB) $(GISLIB) $(PROJLIB)
DEPENDENCIES = $(GPROJDEP) $(RASTERDEP) $(GISDEP)

EXTRA_INC = $(PROJINC)

include $(MODULE_TOPDIR)/include/Make/Module.make

default: cmd
