## adaptiveDijkstra ##
Adaptive Dijkstra algorithm for the horizontal free-flight problem. 

## Installing adaptiveDijkstra ##

A compiler that supports C++11 is needed to build adaptiveDijkstra. Development of the code is performed using version 7.3.0. All dependencies are bundled with adaptiveDijkstra. You can download and compile the latest code from github as follows:

```
git clone --recursive https://github.com/vanhoan310/horizontal_free_flight.git
cd horizontal_free_flight
g++ -std=c++11 -o adaptiveDijkstra -I"include" main.cpp
```

## Running adaptiveDijkstra ##

To generate a sample dataset, e.g filename.txt
./adaptiveDijkstra -gendata filename.txt

To run Adaptive Dijkstra algorithm on data filename.txt
./adaptiveDijkstra filename.txt
