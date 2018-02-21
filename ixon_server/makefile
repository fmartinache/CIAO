all:
	g++ -o ixon_server ixon_server.cpp ImageCreate.c -landor -lm -lpthread

server: server_20171017.o ImageCreate.o
	g++ -o ixon_server ixon_server.o ImageCreate.o -landor -lcfitsio -lpthread -lm

ixon_server.o: ixon_server.cpp
	g++ -c -Wall ixon_server.cpp

ImageCreate.o: ImageCreate.c
	gcc -c ImageCreate.c

clean:
	rm -rf *.o
