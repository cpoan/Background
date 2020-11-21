all: main
main: main.cc
	g++ -g -o 3region4 main.cc -lEG `root-config --cflags` `root-config --libs` -lTreePlayer
hsiaomain: combinedFit.C
	g++ -g -o hsiaomain combinedFit.C -lEG `root-config --cflags` `root-config --libs` -lTreePlayer

clean:
	rm ./*~ ./*.o ./main
