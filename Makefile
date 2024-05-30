CC=gcc

CFLAGS=-Wall  

LIBS= -lgsl

OBJ= main.o bcc_structure.o fcc_structure.o hcp_structure.o Cell_list.o Neigh_list.o Initialize_atom_velocity.o rescale_atom_velocity.o force_calculation.o velocity_verlet.o Euler_forward.o Leap_frog.o Potential.o Total_energy.o stress_calculation.o temperature.o strain.o 
LDFLAGS=-lm 


# phony targets will always be remade, so a file named "clean"

# won't prevent the clean target from running

.PHONY: all clean test


all: main



main: $(OBJ)
	$(CC) -o main $(OBJ) $(CFLAGS) $(LDFLAGS)



clean:
	rm -f main *.o



test: main
	./main

