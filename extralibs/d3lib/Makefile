OUT=../libdftd3.a
OBJS1=f3.o dftd3.o copyc6.o hbsimple.o
OBJS = $(OBJS1) $(OBJS2)

include ../make.in


.PHONY: all
.PHONY: clean
all: $(OUT)


#--------------------------------------------------------------------------
# implizite Regel zur Erzeugung von *.o aus *.F ausschalten
%.o: %.F

# aus *.F mache ein *.f
%.f: %.F
	@echo "making $@ from $<"
	$(CC) $(PREFLAG) $(PREOPTS) $< -o $@

# aus *.f mache ein *.o
%.o: %.f
	@echo "making $@ from $<"
	$(FC) $(FFLAGS) -c $< -o $@

%.o: %.f90
	@echo "making $@ from $<"
	$(FC) $(FFLAGS) -c $< -o $@

# aus *.c mache ein *.o
%.o: %.c
	@echo "making $@ from $<"
	$(CC) $(CCFLAGS) -c $< -o $@

# linken
$(OUT): $(OBJS) 
	ar rcs $(OUT) $(OBJS)

#aufraeumen
clean:
	rm -f *.o $(OUT) *.out

