## Multi-level A* algorithm ##
Multi-level A* algorithm for the horizontal free-flight problem. 

## Installing the package ##

A compiler that supports C++11 is needed to build the package. Development of the code is performed using version 7.3.0. The current version dependencies is [LEMON](https://lemon.cs.elte.hu/trac/lemon). You can download and compile the latest code from github as follows:

```
git clone --recursive https://github.com/vanhoan310/horizontal_free_flight.git
cd horizontal_free_flight
make
```

## Running adaptiveDijkstra ##

To generate a grid graph , e.g., testdata.txt
```
./mastar -gendata testdata.txt
```

To run the algorithm on data testdata.txt
```
./mastar testdata.txt
```
