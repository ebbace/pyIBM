FORT = f77
#FORT = f90
#FLAG = -g -dbl
FLAG = -g

CG2 = d9jr3sp.o
UTL = readfact.o readin.o attn.o option.o
HW  = cgu3hw.o u3mult.o klw.o homgauss.o orthon.o 
SHW = cgu3hwdrv.o $(HW) $(CG2) $(UTL)
SU2 = cgu3drv.o cgu3.o chkwght.o genwght.o $(HW) $(CG2) $(UTL)
SO3 = cgu3o3drv.o cgu3o3.o genL.o dmat.o ovrlap.o kmat.o eispack.o dwr3.o \
      $(HW) $(CG2) $(UTL) readtab.o

KMX = kmatdrv.o kmatsp.o kmat.o dlogf2.o dmat.o eispack.o $(UTL) readtab.o
OVR = ovrldrv.o ovrlsp.o ovrlap.o dlogf2.o kmat.o dmat.o eispack.o $(CG2) $(UTL)
OV1 = ovrl1drv.o ovrlsp.o ovrlap.o dlogf2.o kmat.o dmat.o eispack.o $(CG2) $(UTL)

# -----------------------------------------------------------------------------
#
# SU(2) Clebsch-Gordan coefficients and recoupling coefficients.
#
# -----------------------------------------------------------------------------
dwr3.o      :  dwr3.f
	$(FORT) -c $(FLAG) dwr3.f

d9jr3sp.o   :  d9jr3sp.f
	$(FORT) -c $(FLAG) d9jr3sp.f

# -----------------------------------------------------------------------------
#
# Utilities subprograms.
#
# -----------------------------------------------------------------------------
attn.o      :  attn.f
	$(FORT) -c $(FLAG) attn.f

option.o    :  option.f
	$(FORT) -c $(FLAG) option.f

readfact.o  :  readfact.f
	$(FORT) -c $(FLAG) readfact.f

readin.o    :  readin.f
	$(FORT) -c $(FLAG) readin.f

# -----------------------------------------------------------------------------
#
# SU(3) Clebsch-Gordan coefficients for the highest weight coupled states.
#
# -----------------------------------------------------------------------------
cgu3hwdrv   : $(SHW)
	$(FORT) -o cgu3hwdrv $(FLAG) $(SHW)

cgu3hwdrv.o :  cgu3hwdrv.f
	$(FORT) -c $(FLAG) cgu3hwdrv.f

cgu3hw.o    :  cgu3hw.f
	$(FORT) -c $(FLAG) cgu3hw.f

u3mult.o    :  u3mult.f
	$(FORT) -c $(FLAG) u3mult.f

klw.o       :  klw.f
	$(FORT) -c $(FLAG) klw.f

homgauss.o  :  homgauss.f
	$(FORT) -c $(FLAG) homgauss.f

orthon.o    :  orthon.f
	$(FORT) -c $(FLAG) orthon.f

# -----------------------------------------------------------------------------
#
# SU(3) Clebsch-Gordan coefficients in SU(2) bases.
#
# -----------------------------------------------------------------------------
cgu3drv     : $(SU2) 
	$(FORT) -o cgu3drv $(FLAG) $(SU2)

cgu3drv.o   :  cgu3drv.f
	$(FORT) -c $(FLAG) cgu3drv.f

chkwght.o   :  chkwght.f
	$(FORT) -c $(FLAG) chkwght.f

genwght.o   :  genwght.f
	$(FORT) -c $(FLAG) genwght.f

cgu3.o      :  cgu3.f
	$(FORT) -c $(FLAG) cgu3.f

# -----------------------------------------------------------------------------
#
# SU(3) Clebsch-Gordan coefficients in SO(3) bases.
#
# -----------------------------------------------------------------------------
cgu3o3drv   : $(SO3) Makefile
	$(FORT) -o cgu3o3drv $(FLAG) $(SO3)

cgu3o3drv2  : $(SO3) Makefile
	$(FORT) -o cgu3o3drv2 $(FLAG) $(SO3)

cgu3o3drv.o : cgu3o3drv.f
	$(FORT) -c $(FLAG) cgu3o3drv.f

dmat.o      : dmat.f
	$(FORT) -c $(FLAG) dmat.f

ovrlap.o    : ovrlap.f
	$(FORT) -c $(FLAG) ovrlap.f

kmat.o      : kmat.f
	$(FORT) -c $(FLAG) kmat.f

cgu3o3.o    : cgu3o3.f
	$(FORT) -c $(FLAG) cgu3o3.f

eispack.o   : eispack.f
	$(FORT) -c $(FLAG) eispack.f

genL.o      :  genL.f
	$(FORT) -c $(FLAG) genL.f

# -----------------------------------------------------------------------------
#
# K-matrix for a given SU(3) irrep.
#
# -----------------------------------------------------------------------------
kmatdrv     :  $(KMX) Makefile
	$(FORT) -o kmatdrv $(FLAG) $(KMX)

kmatdrv.o   :  kmatdrv.f
	$(FORT) -c $(FLAG) kmatdrv.f

kmatsp.o    :  kmatsp.f
	$(FORT) -c $(FLAG) kmatsp.f

# -----------------------------------------------------------------------------
#
# Overlap between SU(2) and SO(3) bases of SU(3).
#
# -----------------------------------------------------------------------------
ovrldrv     :  $(OVR) Makefile
	$(FORT) -o ovrldrv $(FLAG) $(OVR)

ovrl1drv    :  $(OV1) Makefile
	$(FORT) -o ovrl1drv $(FLAG) $(OV1)

ovrldrv.o   :  ovrldrv.f
	$(FORT) -c $(FLAG) ovrldrv.f

ovrl1drv.o   :  ovrl1drv.f
	$(FORT) -c $(FLAG) ovrl1drv.f

ovrlsp.o    :  ovrlsp.f
	$(FORT) -c $(FLAG) ovrlsp.f

# -----------------------------------------------------------------------------
#
# Generate all relevant programs
#
# -----------------------------------------------------------------------------
all         :
	make cgu3drv
	make cgu3o3drv
# -----------------------------------------------------------------------------
#
# Clean up the object files
#
# -----------------------------------------------------------------------------
clean       :
	rm *.o
cleanall    :
	rm *.o cgu3drv cgu3o3drv
