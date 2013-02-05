CC=
CFLAGS= 
F90=gfortran
#F90FLAGS=-g -fbounds-check -cpp
F90FLAGS=-g -cpp -O3
PARAFLAGS=-fopenmp

SOURCES= mainvar.f90 \
gen_lib.f90 \
find_defect.f90 \
find_Hbond.f90 \
rdf.f90 \
main.f90 


OBJ=$(addsuffix .o, $(basename $(SOURCES)))


.SUFFIXES :.c .f90

.f90.o:
	$(F90) $(F90FLAGS) $(PARAFLAGS) -c  $< 

.c.o:
	$(CC) $(CFLAGS) -c $<



H3O.x: $(OBJ)
	@echo "Building new Hbond ... "
	@echo "Current objects: $(OBJ)"
	$(F90) $(F90FLAGS) $(PARAFLAGS) $(OBJ) -o $@


clean:	
	@echo "Cleaning Directory ... "
	rm -f $(OBJ) 
 
veryclean: 
	@echo "Cleaning Directory AND Executables ..."
	rm -f $(OBJ) Hbond *.mod 
