// To compile codes
g++ -std=c++11 -o adaptiveDijkstra -I"include" main.cpp

// To generate a sample dataset, e.g filename.txt
./adaptiveDijkstra -gendata filename.txt

// To run Adaptive Dijkstra algorithm on data filename.txt
./adaptiveDijkstra filename.txt
