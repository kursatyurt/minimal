CXX = g++-13 -g -Wall -Wfatal-errors -O3 -fopenmp

all:
	@make kernel
	@make fmm

kernel: kernel.cxx
	$(CXX) $? -o $@
	./kernel 10
	./kernel 20
	./kernel 30
	./kernel 40
	./kernel 50
	./kernel 60

fmm: fmm.cxx
	$(CXX) $? -o $@ -DEXAFMM_EAGER
	./fmm
	$(CXX) $? -o $@ -DEXAFMM_LAZY
	./fmm

clean:
	$(RM) ./*.o ./kernel ./fmm
