all: main.cc
	g++ -g -o unified5 main.cc -lEG `root-config --cflags` `root-config --libs` -lTreePlayer
hsiaomain: combinedFit.C
	g++ -g -o hsiaomain combinedFit.C -lEG `root-config --cflags` `root-config --libs` -lTreePlayer

clean:
	rm ./*~ ./*.o ./main
