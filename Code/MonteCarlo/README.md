To run code, first compile the Jenkins-Traub root finding algorithm with 

g++ -c rpoly_ak1.cpp -O3

Then, compile the running model code,

gcc -c running_model.c -O3

Then link all the object files

g++ -o run.out rpoly_ak1.o running_model.o -lm

Finally, run (this may take a few minutes).

./run.out.

If you wish to run this in the background, use nohup

nohup ./run.out &