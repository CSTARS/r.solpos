MODULE_TOPDIR = ../..

PGM = r.solpos

LIBES = $(GPROJLIB) $(RASTERLIB) $(GISLIB)
DEPENDENCIES = $(GPROJDEP) $(RASTERDEP) $(GISDEP)
EXTRA_INC = $(PROJINC) $(GDALCFLAGS)

include $(MODULE_TOPDIR)/include/Make/Module.make

default: cmd
