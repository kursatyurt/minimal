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
	./kernel 70
	./kernel 80
	./kernel 90
	./kernel 100
	./kernel 110
	./kernel 120
	./kernel 130

fmm: fmm.cxx
	$(CXX) $? -o $@ -DEXAFMM_EAGER
	./fmm

clean:
	$(RM) ./*.o ./kernel ./fmm
