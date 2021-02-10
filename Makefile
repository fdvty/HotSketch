CPFLAGS = -g -rdynamic -mavx2 -mbmi -mavx512bw -mavx512dq --std=c++17 -O3

hotsketch: main.cpp HotSketch.h
	g++ -mavx2 main.cpp -o hotsketch $(CPFLAGS)

clean: 
	rm hotsketch

