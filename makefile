FLAGS=-std=c++11 -O3 -w

main: main.cpp PolyDecomp.o convex_hull.o ArrayX.o
	mpic++ $(FLAGS) $^ -o $@
PolyDecomp.o: PolyDecomp.cpp PolyDecomp.hpp
	mpic++ $(FLAGS) PolyDecomp.cpp PolyDecomp.hpp -c
ArrayX.o: ArrayX.cpp ArrayX.hpp
	g++ $(FLAGS) $^ -c
convex_hull.o: convex_hull.cpp convex_hull.hpp
	g++ $(FLAGS) convex_hull.cpp convex_hull.hpp -c
run:
	mpirun -np 6 -hostfile hostfile main 100 100
clean:
	rm *.o *.gch main
