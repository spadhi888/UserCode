OBJS1 = suspect2_call.o suspect2.o twoloophiggs.o bsg.o
OBJS2 = suspect2.o twoloophiggs.o bsg.f sdecay.o
OBJS3 = susylha.o hgaga.o hdecay.o

FC = g77

.f.o: 
	$(FC) -c -finit-local-zero $*.f


all: suspect2 sdecay

suspect2: $(OBJS1)
	  $(FC) $(OBJS1) -o suspect


sdecay: $(OBJS2) $(OBJS3)
	$(FC) $(OBJS2) $(OBJS3) -o sdecay
