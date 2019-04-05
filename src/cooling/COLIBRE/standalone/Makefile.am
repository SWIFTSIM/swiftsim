CC = gcc

OPTIMIZE = -g 

EXEC   = cool

OBJS = cooling.o cooling_tables.o interpolate.o cooling_rates.o cooling_function.o main.o

INCL = cooling_struct.h cooling_tables.h interpolate.h cooling_rates.h physical_constants_cgs.h cooling_function.h 

HDF5LIB = -lhdf5
LIBS    = -lm $(HDF5LIB) -lgslcblas -lgsl  

$(EXEC): $(OBJS)
	$(CC) $(OPTIMIZE) $(OBJS) $(LIBS) -o $(EXEC)
$(OBJS): $(INCL) 

clean:
	rm -f $(OBJS) $(EXEC) *.i *.s *.o
