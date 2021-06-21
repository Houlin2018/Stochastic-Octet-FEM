octet: main.o functions_Geometry.o functions_Kernel.o functions_Post-processing.o functions_Random_vec.o functions_Step_control.o
	g++ main.o functions_Geometry.o functions_Kernel.o functions_Post-processing.o functions_Random_vec.o functions_Step_control.o -o octet
	
main.o: main.cpp
	g++ -c main.cpp -std=c++11 -I /home/wlm759/Octet_lattice/Eigen3.3.5

functions_Geometry.o: functions_Geometry.cpp functions_Geometry.h
	g++ -c functions_Geometry.cpp -std=c++11
	
functions_Kernel.o: functions_Kernel.cpp functions_Kernel.h
	g++ -c functions_Kernel.cpp -std=c++11 -I /home/wlm759/Octet_lattice/Eigen3.3.5
	
functions_Post-processing.o: functions_Post-processing.cpp functions_Post-processing.h
	g++ -c functions_Post-processing.cpp -std=c++11 -I /home/wlm759/Octet_lattice/Eigen3.3.5
	
functions_Random_vec.o: functions_Random_vec.cpp functions_Random_vec.h
	g++ -c functions_Random_vec.cpp -std=c++11 -I /home/wlm759/Octet_lattice/Eigen3.3.5
	
functions_Step_control.o: functions_Step_control.cpp functions_Step_control.h
	g++ -c functions_Step_control.cpp -std=c++11 -I /home/wlm759/Octet_lattice/Eigen3.3.5
	
clean:
	rm *.o octet