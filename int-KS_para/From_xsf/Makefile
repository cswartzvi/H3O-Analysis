CC=gcc
CFLAGS= 
F90=ftn -mp
#F90FLAGS=-g -traceback -Mbounds
F90FLAGS=
PARAFLAGS=

SOURCES= ks_lib.f90, main.f90


OBJ=$(addsuffix .o, $(basename $(SOURCES)))


.SUFFIXES :.c .f90

.f90.o:
	$(F90) $(F90FLAGS) $(PARAFLAGS) -c  $< 

.c.o:
	$(CC) $(CFLAGS) -c $<



H3O_KS: $(OBJ)
	@echo "Building new KS ... "
	@echo "Current objects: $(OBJ)"
	$(F90) $(F90FLAGS) $(PARAFLAGS) $(OBJ) -o $@


clean:	
	@echo "Cleaning Directory ... "
	rm -f $(OBJ) 
 
veryclean: 
	@echo "Cleaning Directory AND Executables ..."
	rm -f $(OBJ) KS *.mod 
