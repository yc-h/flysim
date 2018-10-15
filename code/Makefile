

flysim.out: Info.o receptor.o StiExt.o neuron_affi.o neuron_gp4.o neuronGnl.o neuronGnl2.o EqsLIF_gp4.o EqsGnl.o EqsGnl1_1.o EqsGnl1_1RK4.o EqsGnl1_1IpvEul.o EqsGnl1_1Eul.o EqsGnl2.o EqsGnl2_1.o parser_affi.o parser.o rbmq.o NetSim.o NetSimGnl.o NetSimGnl1_1.o NetSimGnl1_1_info.o NetSimGnl1_1_BulPp.o NetSimGnl2.o NetSimGnl2_info.o NetSimGnl2_BulPp.o NetSimGnl2_1.o NetSimGnl2_1_info.o NetSimGnl2_1_BulPp.o NetSim_gp4.o NetSimLEGWorm.o main.o

	g++ Info.o receptor.o StiExt.o neuron_affi.o neuron_gp4.o neuronGnl.o neuronGnl2.o EqsLIF_gp4.o EqsGnl.o EqsGnl1_1.o EqsGnl1_1RK4.o EqsGnl1_1IpvEul.o EqsGnl1_1Eul.o EqsGnl2.o EqsGnl2_1.o parser_affi.o parser.o rbmq.o NetSim.o NetSimGnl.o NetSimGnl1_1.o NetSimGnl1_1_info.o NetSimGnl1_1_BulPp.o NetSimGnl2.o NetSimGnl2_info.o NetSimGnl2_BulPp.o NetSimGnl2_1.o NetSimGnl2_1_info.o NetSimGnl2_1_BulPp.o NetSim_gp4.o NetSimLEGWorm.o main.o -o flysim.out -Wall -std=c++11 -pthread -O3

main.o: main.cpp NetSim.h
	g++ -c main.cpp -Wall -std=c++11 -pthread -O3

NetSimLEGWorm.o: NetSimLEGWorm.cpp NetSim.h
	g++ -c NetSimLEGWorm.cpp -Wall -std=c++11 -pthread -O3

NetSim_gp4.o: NetSim_gp4.cpp NetSim.h
	g++ -c NetSim_gp4.cpp -Wall -std=c++11 -pthread -O3

NetSimGnl2_1_BulPp.o: NetSimGnl2_1_BulPp.cpp NetSim.h
	g++ -c NetSimGnl2_1_BulPp.cpp -Wall -std=c++11 -pthread -O3

NetSimGnl2_1_info.o: NetSimGnl2_1_info.cpp NetSim.h
	g++ -c NetSimGnl2_1_info.cpp -Wall -std=c++11 -pthread -O3

NetSimGnl2_1.o: NetSimGnl2_1.cpp NetSim.h
	g++ -c NetSimGnl2_1.cpp -Wall -std=c++11 -pthread -O3

NetSimGnl2_BulPp.o: NetSimGnl2_BulPp.cpp NetSim.h
	g++ -c NetSimGnl2_BulPp.cpp -Wall -std=c++11 -pthread -O3

NetSimGnl2_info.o: NetSimGnl2_info.cpp NetSim.h
	g++ -c NetSimGnl2_info.cpp -Wall -std=c++11 -pthread -O3

NetSimGnl2.o: NetSimGnl2.cpp NetSim.h
	g++ -c NetSimGnl2.cpp -Wall -std=c++11 -pthread -O3

NetSimGnl1_1_BulPp.o: NetSimGnl1_1_BulPp.cpp NetSim.h
	g++ -c NetSimGnl1_1_BulPp.cpp -Wall -std=c++11 -pthread -O3

NetSimGnl1_1_info.o: NetSimGnl1_1_info.cpp NetSim.h
	g++ -c NetSimGnl1_1_info.cpp -Wall -std=c++11 -pthread -O3

NetSimGnl1_1.o: NetSimGnl1_1.cpp NetSim.h
	g++ -c NetSimGnl1_1.cpp -Wall -std=c++11 -pthread -O3

NetSimGnl.o: NetSimGnl.cpp NetSim.h
	g++ -c NetSimGnl.cpp -Wall -std=c++11 -pthread -O3

NetSim.o: NetSim.cpp NetSim.h
	g++ -c NetSim.cpp -Wall -std=c++11 -pthread -O3

rbmq.o: rbmq.cpp NetSim.h
	g++ -c rbmq.cpp -Wall -std=c++11 -pthread -O3

parser.o: parser.cpp parser.h
	g++ -c parser.cpp -Wall -std=c++11 -pthread -O3

parser_affi.o: parser_affi.cpp parser.h
	g++ -c parser_affi.cpp -Wall -std=c++11 -pthread -O3

EqsGnl2_1.o: EqsGnl2_1.cpp Eqs.h
	g++ -c EqsGnl2_1.cpp -Wall -std=c++11 -pthread -O3

EqsGnl2.o: EqsGnl2.cpp Eqs.h
	g++ -c EqsGnl2.cpp -Wall -std=c++11 -pthread -O3

EqsGnl1_1RK4.o: EqsGnl1_1RK4.cpp Eqs.h
	g++ -c EqsGnl1_1RK4.cpp -Wall -std=c++11 -pthread -O3

EqsGnl1_1IpvEul.o: EqsGnl1_1IpvEul.cpp Eqs.h
	g++ -c EqsGnl1_1IpvEul.cpp -Wall -std=c++11 -pthread -O3

EqsGnl1_1Eul.o: EqsGnl1_1Eul.cpp Eqs.h
	g++ -c EqsGnl1_1Eul.cpp -Wall -std=c++11 -pthread -O3

EqsGnl1_1.o: EqsGnl1_1.cpp Eqs.h
	g++ -c EqsGnl1_1.cpp -Wall -std=c++11 -pthread -O3

EqsGnl.o: EqsGnl.cpp Eqs.h
	g++ -c EqsGnl.cpp -Wall -std=c++11 -pthread -O3

EqsLIF_gp4.o: EqsLIF_gp4.cpp Eqs.h
	g++ -c EqsLIF_gp4.cpp -Wall -std=c++11 -pthread -O3

neuronGnl2.o: neuronGnl2.cpp neuron.h
	g++ -c neuronGnl2.cpp -Wall -std=c++11 -pthread -O3

neuronGnl.o: neuronGnl.cpp neuron.h
	g++ -c neuronGnl.cpp -Wall -std=c++11 -pthread -O3

neuron_gp4.o: neuron_gp4.cpp neuron.h
	g++ -c neuron_gp4.cpp -Wall -std=c++11 -pthread -O3

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




