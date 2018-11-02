

flysim.out: Info.o receptor.o StiExt.o neuron_affi.o neuronGnl2.o EqsGnl2.o parser_affi.o parser.o NetSim.o NetSimGnl2.o NetSimGnl2_info.o NetSimGnl2_BulPp.o main.o

	g++ Info.o receptor.o StiExt.o neuron_affi.o neuronGnl2.o EqsGnl2.o parser_affi.o parser.o NetSim.o NetSimGnl2.o NetSimGnl2_info.o NetSimGnl2_BulPp.o main.o -o flysim.out -Wall -std=c++11 -pthread -O3

main.o: main.cpp NetSim.h
	g++ -c main.cpp -Wall -std=c++11 -pthread -O3

NetSimGnl2_BulPp.o: NetSimGnl2_BulPp.cpp NetSim.h
	g++ -c NetSimGnl2_BulPp.cpp -Wall -std=c++11 -pthread -O3

NetSimGnl2_info.o: NetSimGnl2_info.cpp NetSim.h
	g++ -c NetSimGnl2_info.cpp -Wall -std=c++11 -pthread -O3

NetSimGnl2.o: NetSimGnl2.cpp NetSim.h
	g++ -c NetSimGnl2.cpp -Wall -std=c++11 -pthread -O3

NetSim.o: NetSim.cpp NetSim.h
	g++ -c NetSim.cpp -Wall -std=c++11 -pthread -O3

parser.o: parser.cpp parser.h
	g++ -c parser.cpp -Wall -std=c++11 -pthread -O3

parser_affi.o: parser_affi.cpp parser.h
	g++ -c parser_affi.cpp -Wall -std=c++11 -pthread -O3

EqsGnl2.o: EqsGnl2.cpp Eqs.h
	g++ -c EqsGnl2.cpp -Wall -std=c++11 -pthread -O3

neuronGnl2.o: neuronGnl2.cpp neuron.h
	g++ -c neuronGnl2.cpp -Wall -std=c++11 -pthread -O3

neuron_affi.o: neuron_affi.cpp neuron_affi.h
	g++ -c neuron_affi.cpp -Wall -std=c++11 -pthread -O3

StiExt.o: StiExt.cpp StiExt.h
	g++ -c StiExt.cpp -Wall -std=c++11 -pthread -O3

receptor.o: receptor.cpp receptor.h
	g++ -c receptor.cpp -Wall -std=c++11 -pthread -O3

Info.o: Info.cpp Info.h
	g++ -c Info.cpp -Wall -std=c++11 -pthread -O3

clean:
	rm *.o *.out *~




