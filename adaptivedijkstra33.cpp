#include "adaptivedijkstra33.h"
#include <iomanip>
#include <math.h>
#include <vector>
#include <assert.h>
#include <lemon/dim2.h>
#include <lemon/list_graph.h>
#include <lemon/adaptors.h>  
#include <lemon/dijkstra.h>
#include <unordered_map>
#include "filterastarv19.h" 
#include <fstream>
#include <cstdlib>
#include <sstream> 
#include <string>  
using namespace lemon;
using namespace std;
typedef ListDigraph Graph;
typedef Graph::Node Node;
typedef Graph::Arc Arc;
typedef Graph::ArcMap<double> LengthMap;
typedef Graph::ArcIt ArcIt;
typedef Graph::NodeIt NodeIt;
AdaptDijkstra::AdaptDijkstra (Graph & g, LengthMap& length, int numberofStage, Graph::NodeMap<Point> & coordinate, int numxbase, int numybase):
	_graph(g), _length(length), _stage(1),_distSource(g), _distSink(g),potential(g,_INFINITY),
	_filter(g,false),_height(g), _maxStage(numberofStage-1),_neighbor(g),_parent(g),_nodeNumber(countNodes(_graph)),
	_numxbase(numxbase), _numybase(numybase),_runtime(0),_coordinate(coordinate)
	/* numberofStage: number of stages needed to process the algorithm.
	stage = 1 tuong ung voi do thi thua nhat, stage = 2 do thi day hon 1 cap, ..,
	stage = numberofStage  tuong ung voi do thi min nhat G. */
	{}
	// Compute power of integer: base^exp //TODO
inline double AdaptDijkstra::calculateHeuristic (ListDigraph::NodeMap<double> &dist,  ListDigraph::NodeMap<vector<Node>> & _parent, ListDigraph::Node & v)
{
	if(_parent[v].at(0) == lemon::INVALID) return min(dist[v],2000.0);	
	//cout << _parent[v].size() << "\t";
	double heuristicValue = dist[v];
	for (int i = 0; i< _parent[v].size(); ++i) 
	{  
		Node w = _parent[v].at(i);
		if (dist[w] == _INFINITY) dist[w] = calculateHeuristic(dist,_parent,w);
		if (dist[w] < heuristicValue )
			heuristicValue = dist[w];
	}
	dist[v] = heuristicValue; //TODO co nen gan luon gia tri cho thang nay ko ???
	return heuristicValue; 
}    

void AdaptDijkstra::findNeighborNextStage(int currentStage)
{
	const int k = _maxStage - currentStage;
	int nodeNum, nodeNumleft, nodeNumright;
	Node v;
	// Clear all old neighbors
	for (int i = 0; i < _numxbase; i=i+1){ //ipow(2,k-1)
		for (int j = 0; j < _numybase; j=j+1){
		_neighbor[_graph.nodeFromId(i + j*_numxbase)].clear();
	}}
	// Update neighbors for next stage without boundary.//TODO 
	int step_size = ipow(2,k-1);
	int big_step = step_size*2;
	for (int i = big_step; i < _numxbase-1; i=i+big_step){
		for (int j = big_step; j < _numybase-1; j=j+big_step){
		nodeNum     = i + j*_numxbase;
		nodeNumleft     = nodeNum - step_size;
		nodeNumright     = nodeNum+ step_size;
		v = _graph.nodeFromId(nodeNum);
		
// 		_neighbor[v].push_back(v);
		_neighbor[_graph.nodeFromId(nodeNumleft)].push_back(v);
		_neighbor[_graph.nodeFromId(nodeNumright)].push_back(v);
		
		_neighbor[_graph.nodeFromId(nodeNum+_numxbase*step_size)].push_back(v);
		_neighbor[_graph.nodeFromId(nodeNumright+_numxbase*step_size)].push_back(v);
		_neighbor[_graph.nodeFromId(nodeNumleft+_numxbase*step_size)].push_back(v);
		
		_neighbor[_graph.nodeFromId(nodeNum-_numxbase*step_size)].push_back(v);
		_neighbor[_graph.nodeFromId(nodeNumright-_numxbase*step_size)].push_back(v);
		_neighbor[_graph.nodeFromId(nodeNumleft-_numxbase*step_size)].push_back(v);
		
	}}

	// i = 0
		for (int j = ipow(2,k); j < _numybase-1; j=j+ipow(2,k)){
		nodeNum     =  j*_numxbase;
		nodeNumright     = nodeNum+ ipow(2,k-1);
		v = _graph.nodeFromId(nodeNum);
		
// 		_neighbor[v].push_back(v);
		//if (nodeNumright < _nodeNumber )
		_neighbor[_graph.nodeFromId(nodeNumright)].push_back(v);
		//if (nodeNum+ _numxbase < _nodeNumber )
		_neighbor[_graph.nodeFromId(nodeNum+_numxbase*step_size)].push_back(v);
		//if (nodeNumright+ _numxbase < _nodeNumber)
		_neighbor[_graph.nodeFromId(nodeNumright+_numxbase*step_size)].push_back(v);
		//if (nodeNum- _numxbase >= 0)
			_neighbor[_graph.nodeFromId(nodeNum-_numxbase*step_size)].push_back(v);
		//if (nodeNumright- _numxbase >= 0)
			_neighbor[_graph.nodeFromId(nodeNumright-_numxbase*step_size)].push_back(v);
		}
	// i = _numxbase - 1
	for (int j = ipow(2,k); j < _numybase-1; j=j+ipow(2,k)){
		nodeNum     = _numxbase-1 + j*_numxbase;
		nodeNumleft     = nodeNum - ipow(2,k-1);
		v = _graph.nodeFromId(nodeNum);
		
// 		_neighbor[v].push_back(v);
		//if (nodeNumleft >= 0 )
		_neighbor[_graph.nodeFromId(nodeNumleft)].push_back(v);
		//if (nodeNum+ _numxbase < _nodeNumber )
		_neighbor[_graph.nodeFromId(nodeNum+_numxbase*step_size)].push_back(v);
		//if (nodeNumleft+ _numxbase < _nodeNumber)
		_neighbor[_graph.nodeFromId(nodeNumleft+_numxbase*step_size)].push_back(v);
		//if (nodeNum-_numxbase >=0)
		_neighbor[_graph.nodeFromId(nodeNum-_numxbase*step_size)].push_back(v);
		//if (nodeNumleft - _numxbase >=0)
			_neighbor[_graph.nodeFromId(nodeNumleft-_numxbase*step_size)].push_back(v);       
		}
	// j = 0
	for (int i = ipow(2,k); i < _numxbase-1; i=i+ipow(2,k)){
		nodeNum     = i;
		nodeNumleft     = nodeNum - ipow(2,k-1);
		nodeNumright     = nodeNum+ ipow(2,k-1);
		v = _graph.nodeFromId(nodeNum);
		
// 		_neighbor[v].push_back(v);
		//if (nodeNumleft >= 0 )
		_neighbor[_graph.nodeFromId(nodeNumleft)].push_back(v);
		//if (nodeNumright < _nodeNumber )
		_neighbor[_graph.nodeFromId(nodeNumright)].push_back(v);
		
		//if (nodeNum+ _numxbase < _nodeNumber )
		_neighbor[_graph.nodeFromId(nodeNum+_numxbase*step_size)].push_back(v);
		//if (nodeNumright+ _numxbase < _nodeNumber)
		_neighbor[_graph.nodeFromId(nodeNumright+_numxbase*step_size)].push_back(v);
		//if (nodeNumleft+ _numxbase < _nodeNumber)
		_neighbor[_graph.nodeFromId(nodeNumleft+_numxbase*step_size)].push_back(v);
				
	}
	// j = _numxbase-1
	for (int i = ipow(2,k); i < _numxbase-1; i=i+ipow(2,k)){
		nodeNum     = i + (_numxbase-1)*_numxbase;
		nodeNumleft     = nodeNum - ipow(2,k-1);
		nodeNumright     = nodeNum+ ipow(2,k-1);
		v = _graph.nodeFromId(nodeNum);
		
// 		_neighbor[v].push_back(v);
		//if (nodeNumleft >= 0 )
		_neighbor[_graph.nodeFromId(nodeNumleft)].push_back(v);
		//if (nodeNumright < _nodeNumber )
		_neighbor[_graph.nodeFromId(nodeNumright)].push_back(v);
		_neighbor[_graph.nodeFromId(nodeNum-_numxbase*step_size)].push_back(v);
		_neighbor[_graph.nodeFromId(nodeNumright-_numxbase*step_size)].push_back(v);
		_neighbor[_graph.nodeFromId(nodeNumleft-_numxbase*step_size)].push_back(v);
	}
	// i = 0, j = 0
	v = _graph.nodeFromId(0);
// 	_neighbor[v].push_back(v);
	_neighbor[_graph.nodeFromId(step_size)].push_back(v);
	_neighbor[_graph.nodeFromId(_numxbase*step_size)].push_back(v);
	_neighbor[_graph.nodeFromId(step_size+_numxbase*step_size)].push_back(v);
	// i = 0, j = _numxbase
	nodeNum     = (_numxbase-1)*_numxbase;
	v = _graph.nodeFromId(nodeNum);
// 	_neighbor[v].push_back(v);
	_neighbor[_graph.nodeFromId(nodeNum+step_size)].push_back(v);
	_neighbor[_graph.nodeFromId(nodeNum-_numxbase*step_size)].push_back(v);
	_neighbor[_graph.nodeFromId(nodeNum+step_size-_numxbase*step_size)].push_back(v);
	// i = _numxbase, j = 0;
	nodeNum     = _numxbase-1;
	v = _graph.nodeFromId(nodeNum);
// 	_neighbor[v].push_back(v);
	_neighbor[_graph.nodeFromId(nodeNum-step_size)].push_back(v);
	_neighbor[_graph.nodeFromId(nodeNum+_numxbase*step_size)].push_back(v);
	_neighbor[_graph.nodeFromId(nodeNum-step_size+_numxbase*step_size)].push_back(v);
	// i = _numxbase -1, j = i
	nodeNum     = _numxbase-1 + (_numxbase-1)*_numxbase;
	v = _graph.nodeFromId(nodeNum);
// 	_neighbor[v].push_back(v);
	_neighbor[_graph.nodeFromId(nodeNum-step_size)].push_back(v);
	_neighbor[_graph.nodeFromId(nodeNum-_numxbase*step_size)].push_back(v);
	_neighbor[_graph.nodeFromId(nodeNum-step_size-_numxbase*step_size)].push_back(v);
    
}


void AdaptDijkstra::findNeighborNextStageSmall(int currentStage)
{
	const int k = _maxStage - currentStage;
	int nodeNum, nodeNumleft, nodeNumright;
	Node v;
	// Clear all old neighbors
	for (Graph::NodeIt v(_graph); v != INVALID; ++v){
		_neighbor[v].clear();		
	}
// 	for (int i = 0; i < _numxbase; i=i+1){ //ipow(2,k-1)
// 		for (int j = 0; j < _numybase; j=j+1){
// 		_neighbor[_graph.nodeFromId(i + j*_numxbase)].clear();
// 	}}
	// Update neighbors for next stage without boundary.//TODO 
	int step_size = ipow(2,k-1);
	int big_step = step_size*2;
	for (int i = big_step; i < _numxbase-1; i=i+big_step){
		for (int j = big_step; j < _numybase-1; j=j+big_step){
		nodeNum     = i + j*_numxbase;
		nodeNumleft     = nodeNum - step_size;
		nodeNumright     = nodeNum+ step_size;
		v = _graph.nodeFromId(nodeNum);
		
		
		_neighbor[_graph.nodeFromId(nodeNumleft)].push_back(v);
		_neighbor[_graph.nodeFromId(nodeNumright)].push_back(v);
		
		_neighbor[_graph.nodeFromId(nodeNum+_numxbase*step_size)].push_back(v);
		_neighbor[_graph.nodeFromId(nodeNumright+_numxbase*step_size)].push_back(v);
		_neighbor[_graph.nodeFromId(nodeNumleft+_numxbase*step_size)].push_back(v);
		
		_neighbor[_graph.nodeFromId(nodeNum-_numxbase*step_size)].push_back(v);
		_neighbor[_graph.nodeFromId(nodeNumright-_numxbase*step_size)].push_back(v);
		_neighbor[_graph.nodeFromId(nodeNumleft-_numxbase*step_size)].push_back(v);
		
	}}
}

void AdaptDijkstra::constructParentNodes()
{
	for (int i = 0; i < _numxbase; i=i+ipow(2,_maxStage-1)){
		for (int j = 0; j < _numybase; j=j+ipow(2,_maxStage-1)){
		_parent[_graph.nodeFromId(i + j*_numxbase)].push_back(lemon::INVALID);
	}}
	
	for (int stage = _maxStage - 2; stage < _maxStage; ++stage)
	{
		
		/*if(stage < _maxStage - 1) findNeighborNextStageSmall(stage);
		else findNeighborNextStage(stage); */
		findNeighborNextStageSmall(stage);
		for (int i = 0; i < _numxbase; i=i+ipow(2,_maxStage-stage-1)){
		for (int j = 0; j < _numybase; j=j+ipow(2,_maxStage-stage-1)){
			Node v = _graph.nodeFromId(i + j*_numxbase);
//  			_parent[v].push_back(v);  //TODO
			for (int k = 0; k < _neighbor[v].size() ; ++k)
			{
				_parent[v].push_back(_neighbor[v].at(k));
			}
		}}		
	}
}

void AdaptDijkstra::init (Node s, Node t) // Initialize the data
{
	
	// Set _d[v] = _INFINITY for every vertices v of graph G.
	for (Graph::NodeIt v(_graph); v != INVALID; ++v){
		(_distSource)[v] = _INFINITY; (_distSink)[v] = _INFINITY;  		
	}
	_source = s; _sink = t;
// 	for (Graph::NodeIt v(_graph); v != INVALID; ++v){
// 	(_distSource)[v] = calculateDistance(_coordinate,s,v); (_distSink)[v] = calculateDistance(_coordinate,t,v);   }
// // 	
	// Construct height function	
	for (int stage = _maxStage -1 ; stage >= 0; --stage)
	{
		for (int i = 0; i < _numxbase; i=i+ipow(2,_maxStage-stage-1)){
		for (int j = 0; j < _numybase; j=j+ipow(2,_maxStage-stage-1)){
			Node v = _graph.nodeFromId(i + j*_numxbase);
			_height[v] = stage;
		}}	
	}
	
	_height[s] = 0; _height[t] = 0;
	
}

void AdaptDijkstra::start (Node s, Node t)  // Start the alg in the first stage
{
	//- Run Dijkstra on graph H from source = s to every nodes. 2^(_maxStage - 1)
	
	astar2dir(filterNodes( _graph, lemon::lessMap(_height, constMap<ListDigraph::Node>(_stage))) ,_length).path(_shortest_path).
	  dist(_shortest_distance).runastarAdaptive(s,t,_distSource, 0.1,_INFINITY); //TODO 0.1 means 10% of overun dijkstra
	  _adaptFactor = 1.01;

	//- Increase the _stage by 1;
	_stage ++;
}

void AdaptDijkstra::processNextStage(Node s, Node t, double & time, bool savePathNodes, Graph::NodeMap<vector<Node>> & neighborSet)  // Process the next stage of AdaptDijkstra.
{
	/* In fact, it is not necessary to compute pi(v) for every vertex in graph H. We can put _neighbor(G)
	into A* algorithm as an argument to calculate pi(v) for the needed vertex. */
	//for (int i=0; i< 20;++i) if (_filter[_graph.nodeFromId(i)] == true) cout << i << "\t";
		
	if (_stage % 2 == 0)
	{
		// Run A* REVERSE alg with heuristic function pi(V(H)) from node t to every nodes in H;
		
		/*
		for (Graph::NodeIt v(_graph); v != INVALID; ++v)
		{
			if(_filter[v] == true)
				potential[v] = calculateHeuristic(_distSource,_parent,v);//cout << _graph.id(v) << "t"<< temp_dist[v] <<endl;}
		}
		*/
		clock_t t1,t2;t1=clock();
		double distance;
		astar2dir(reverseDigraph(filterNodes( _graph, lemon::lessMap(_height, constMap<ListDigraph::Node>(_stage)))),_length).path(_shortest_path).
		dist(distance).runastar(t,s,neighborSet,_distSink, _distSource,potential, _shortest_distance, _INFINITY,_parent,_adaptFactor);
		_adaptFactor = _shortest_distance/distance;
		_shortest_distance = distance;
		//for (int i = 0; i < 20; i++) if (_filter[_graph.nodeFromId(i)]==true) cout << i << "\t";
		t2=clock();    float diff ((float)t2-(float)t1);
		time = time + diff/(CLOCKS_PER_SEC);
	}	
	else
	{
		//- Run A* alg with heuristic function pi(V(H)) from node s to every nodes in H;
		double distance;
		clock_t t1,t2;t1=clock();  
		astar2dir(filterNodes( _graph, lemon::lessMap(_height, constMap<ListDigraph::Node>(_stage))),_length).path(_shortest_path).
		dist(distance).runastar(s,t,neighborSet,_distSource, _distSink,potential, _shortest_distance, _INFINITY,_parent, _adaptFactor);
		_adaptFactor = _shortest_distance/distance;
		_shortest_distance = distance;
		//for (int i = 0; i < 20; i++) if (_filter[_graph.nodeFromId(i)]==true) cout << i << "\t";
		t2=clock();    float diff ((float)t2-(float)t1);
		time = time + diff/(CLOCKS_PER_SEC);

	}
	
	// Save the path node
	if(savePathNodes == true)
		{
		for (int i=0; i < _shortest_path.length() ; i++)
		{
			Arc e = _shortest_path.nth(i);
			_old_path_Nodes.push_back(_graph.target(e));
		}
		_old_path_Nodes.push_back(s);
		_old_path_Nodes.push_back(t);
		}
	//- Inscrease the _stage by 1
	_stage ++;
}

void AdaptDijkstra::run(Node s,Node t)
{	

// 	constructParentNodes();
	init(s,t);     // Run the intilization	
	clock_t t1,t2;
	t1=clock();		
 	start(s,t);  // Run the first stage
	t2=clock();    float diff ((float)t2-(float)t1);
	_runtime = _runtime + diff/(CLOCKS_PER_SEC);
  	for(int i=2; i<= _maxStage - 1; ++i )
	{
		findNeighborNextStageSmall(_stage-1);
 		processNextStage(s,t,_runtime,false,_neighbor);  				 
	}
	
// 	processNextStage(s,t,_runtime,false,_neighbor);
// 	processNextStage(s,t,_runtime,false,_parent); 
	
	findNeighborNextStageSmall(_stage-1);
	t1=clock();
	if (_stage % 2 == 0)
	{
		// Run A* REVERSE alg with heuristic function pi(V(H)) from node t to every nodes in H;
		
		/*
		for (Graph::NodeIt v(_graph); v != INVALID; ++v)
		{
			if(_filter[v] == true)
				potential[v] = calculateHeuristic(_distSource,_parent,v);//cout << _graph.id(v) << "t"<< temp_dist[v] <<endl;}
		}
		*/
// 		clock_t t1,t2;t1=clock();
		double distance;
		astar2dir(reverseDigraph(_graph),_length).path(_shortest_path).dist(distance).runastar(t,s,_neighbor,_distSink, _distSource,potential, _shortest_distance, _INFINITY,_parent,_adaptFactor);
		_adaptFactor = _shortest_distance/distance;
		_shortest_distance = distance;
	}	
	else
	{
		double distance;
		astar2dir(_graph,_length).path(_shortest_path).dist(distance).runastar(s,t,_neighbor,_distSource, _distSink,potential, _shortest_distance, _INFINITY,_parent, _adaptFactor);
		_adaptFactor = _shortest_distance/distance;
		_shortest_distance = distance;
	}
    
    //- Inscrease the _stage by 1
	_stage ++;
	t2=clock();    float diff1 (t2-t1);
 	_runtime = _runtime + diff1/(CLOCKS_PER_SEC);
}

/*
vector<Node> AdaptDijkstra::getShortestPath()
{
	Arc e; Node a;
	if (_maxStage % 2 == 1)
	{
		e = _shortest_path.front();
		_path_Nodes.push_back(_graph.source(e));
		while (_shortest_path.empty() == 0)
		{
		e = _shortest_path.front();
		_path_Nodes.push_back(_graph.target(e));
		_shortest_path.eraseFront();
		}
	}
	else {
		e = _shortest_path.back();
		_path_Nodes.push_back(_graph.source(e));
		while (_shortest_path.empty() == 0)
		{
		e = _shortest_path.back();
		_path_Nodes.push_back(_graph.target(e));
		_shortest_path.eraseBack();
		}
	}
	return _path_Nodes;
} */

inline int AdaptDijkstra::ipow(int base, int exp)
{
	int result = 1;
	while (exp)
	{
		if (exp & 1)
		result *= base;
		exp >>= 1;
		base *= base;
	}
	
	return result;
}

double AdaptDijkstra::getShortestDistance()
{
	return _shortest_distance;
}

double AdaptDijkstra::getRunTime()
{
	return _runtime;
}

double AdaptDijkstra::calculateDistance(ListDigraph::NodeMap<lemon::dim2::Point<double>>& coordbase,
	Node currentNode, Node sinkNode) 
{
	const double x1 = coordbase[currentNode].x, y1 = coordbase[currentNode].y;
	const double x2 = coordbase[sinkNode].x,    y2 = coordbase[sinkNode].y;
	return sqrt((x1-x2)*(x1-x2)+(y1-y2)*(y1-y2));
}

void AdaptDijkstra::printShortestPath(std::ofstream& outdata)
{
	outdata << "x =[";  
	outdata << "" << _coordinate[_source].x << "," << _coordinate[_source].y << ";";
	if (_maxStage%2 != 0)
	{
		for (int i= 0; i < _shortest_path.length(); i++)
		{
		Arc e = _shortest_path.nth(i);
		Node a = _graph.source(e);
		outdata << "" << _coordinate[a].x << "," << _coordinate[a].y << ";";
		}
	}
	else{
		for (int i= _shortest_path.length()-1; i >= 0 ; i--)
		{
		Arc e = _shortest_path.nth(i);
		Node a = _graph.target(e);
		outdata << "" << _coordinate[a].x << "," << _coordinate[a].y << ";";
		}
	}
	outdata << "" << _coordinate[_sink].x << "," << _coordinate[_sink].y << "];\n";
	outdata << "hold on " << "\n";
	outdata << "plot(x(:,1),x(:,2),'r','LineWidth',2)" << "\n";
	outdata << " axis([-20000 20000 -10000 10000])" << "\n"; //TODO
	outdata << " text(x(1,1),x(1,2),'from Airport')" << "\n";
	outdata << " text(x(end,1),x(end,2),'to Airport')" << "\n";	  
	outdata << "[x,y] = meshgrid(-20000:400:20000,-10000:400:10000);" << "\n";
	outdata << "zx = -2000;zy = 7000;z2x = 8000;z2y = -6000;" << "\n";
	outdata << "u = ((y-zy).*exp(-((x-zx).^2 + (y-zy).^2)/2000000.0))/2;" << "\n";		
	outdata << "v = ((zx-x).*exp(-((x-zx).^2 + (y-zy).^2)/2000000.0))/2; " << "\n";
	outdata << "quiver(x,y,u,v,1.2) " << "\n";
	outdata << "axis square" << "\n";
	outdata << "print figure.jpg" << "\n";
	outdata << "hold off; " << "\n";
	//outdata << " " << "\n";
	outdata.close();
	return;
}

/**
 * Computes overfly costs corresponding to this airspace
 * @param weight Weight in tonnes (call weightType() first to determine which weight to use)
 * @param distance Distance in km (call distanceType() first to determine which distance to use)
 * @param fromAirportId Lido Id of departure airport
 * @param toAirportId Lido Id of destination airport
 * @return total overflight cost in USD
 */
MakeGraph::MakeGraph (std::string filename, const int numxNodes, const int numyNodes, const double xmin, const double xmax, const double ymin, const double ymax, double radius):
	_filename(filename), _numxNodes(numxNodes), _numyNodes(numyNodes), _xmin(xmin), _xmax(xmax), _ymin(ymin),_ymax(ymax), _radius(radius)
	{}
	
void MakeGraph::getGroundDistance(const double& x1,const double& y1,const double& x2, 
	const double& y2,const int& numx,const int& numy,const double& dx,const double& dy
	,double & arclength, double & arclengthv2, double _xmin, double _ymin)
{
    double minx = (x1<x2)? x1:x2, maxx = (x1<x2)? x2:x1;
    double miny = (y1<y2)? y1:y2, maxy = (y1<y2)? y2:y1; 
    double wCross, wTrack;  
    double vPlane = 850.0;
    double a,b,temp;
    double distance2P = distance2Points(x1,y1,x2,y2);
    double cosAngle = (y2-y1)/distance2P, sinAngle = (x2-x1)/distance2P;
    arclength = 0.0; arclengthv2 = 0.0;     
    if (x1==x2)
    { 
            if (y1==y2) return;
    //         int columnNum = round((x1-_xmin)/dx);
    //         int startRow = round((miny-_ymin)/dy);
            int maxIter = round((maxy-miny)/dy);
            for (int i=0;i< maxIter;++i)
        {
            double wEast  = getWindEast(x1,miny + i*dy+.5*dy);
            double wNorth = getWindNorth(x1,miny + i*dy+.5*dy);
            wCross = -wEast*cosAngle+wNorth*sinAngle;
            wTrack =  wEast*sinAngle+wNorth*cosAngle; 
            temp = sqrt(vPlane*vPlane - wCross*wCross);
            arclength   += dy/(temp+wTrack);
            arclengthv2 += dy/(temp-wTrack);
        }
        arclength = arclength*vPlane;
        arclengthv2 = arclengthv2*vPlane;
        return;
    }
    
    if (y1==y2)
    {
        int rowNum = round(y1/dy);
        int startColumn = round(minx/dx);
        int maxIter = round((maxx-minx)/dx);
        for (int i=0;i<maxIter;++i)
        {		
            double wEast  = getWindEast (y1,minx + i*dx + 0.5*dx);
            double wNorth = getWindNorth(y1,minx + i*dx + 0.5*dx);
            wCross = -wEast*cosAngle+wNorth*sinAngle;
            wTrack =  wEast*sinAngle+wNorth*cosAngle; 
            temp = sqrt(vPlane*vPlane - wCross*wCross);
            arclength += dx/(temp+wTrack);
            arclengthv2 += dx/(temp-wTrack);
        }
        arclength = arclength*vPlane;
        arclengthv2 = arclengthv2*vPlane;
        return;
    }
    
	vector<double> v1, v2, vectorMerged;
	int xIter = round((maxx-minx)/dx);
	int yIter = round((maxy-miny)/dy);
	for (int i=0;i<=xIter;++i) v1.push_back(minx+i*dx);
	
	a = (x2-x1)/(y2-y1); b = x1 - a*y1;
	if (a>0)
	{
		for (int i=0;i<=yIter;++i) v2.push_back(a*(miny+i*dy)+b);
	} 
	else
	{
		for (int i=yIter;i>=0;--i) v2.push_back(a*(miny+i*dy)+b);
	}
	int vectorsize = v1.size() + v2.size();
	vectorMerged.resize(vectorsize);
	merge(v1.begin(), v1.end(), v2.begin(),v2.end(), vectorMerged.begin());
	double px1,px2,py1,py2;
	arclength = 0.0; arclengthv2 = 0.0;
	for (int i=0;i<vectorMerged.size()-1; ++i)
	{
		px1 = vectorMerged.at(i)  ; 
		px2 = vectorMerged.at(i+1); 
		if (px1 != px2)
		{
			py1 = Line(x1,y1,x2,y2,px1);
			py2 = Line(x1,y1,x2,y2,px2);
			double wEast  = getWindEast((px1+px2)/2,(py1+py2)/2);
			double wNorth = getWindNorth((px1+px2)/2,(py1+py2)/2);
	// 	 	cout << wEast << "\t" << wNorth << "\t";
			wCross = -wEast*cosAngle+wNorth*sinAngle;
			wTrack =  wEast*sinAngle+wNorth*cosAngle; 
			temp = sqrt(vPlane*vPlane - wCross*wCross);
	// 		cout << temp<< "\t";
			distance2P = distance2Points(px1,py1,px2,py2);
			arclength = arclength + distance2P/(temp+wTrack);
			arclengthv2 += distance2P/(temp-wTrack);
			std::cout.precision(2);
// 			cout << arclength << "\t";
		}
	}
	arclength = arclength*vPlane;
	arclengthv2 = arclengthv2*vPlane;
}
//v
double MakeGraph::getWindNorth(double x,double y) { return 0.5*((-2000.0-x)*exp(-(pow(x+2000.0,2) + pow(y-7000.0,2))/2000000.0));} 
//u
double MakeGraph::getWindEast(double x,double y)  { return 0.5*((y-7000.0) *exp(-(pow(x+2000.0,2) + pow(y-7000.0,2))/2000000.0));}

std::vector<double> MakeGraph::merge2Sorted( const std::vector<double>& left, const std::vector<double>& right ) {
    std::vector<double> output;
    std::merge(left.begin(), left.end(), right.begin(), right.end(), std::back_inserter(output));
    return output;
}

const double MakeGraph::Line(const double& x1, const double& y1, const double& x2, const double& y2, const double& a)
{
	return (y2-y1)*(a-x1)/(x2-x1)+y1;
}

const double MakeGraph::distance2Points(const double& x1, const double& y1, const double& x2, const double& y2)
{
    return sqrt((x1-x2)*(x1-x2)+(y1-y2)*(y1-y2));
}

void MakeGraph::makeGraphData()
{
	typedef ListDigraph Graph;
	typedef Graph::Node Node;
	typedef Graph::Arc Arc;
	typedef Graph::ArcMap<double> LengthMap;
	typedef Graph::ArcIt ArcIt;
	typedef Graph::NodeIt NodeIt;
	Graph gbase;  
	typedef dim2::Point<double> Point;
	Graph::NodeMap<Point> coordbase(gbase);
	LengthMap lenbase(gbase);
	Graph::NodeMap<int> _height(gbase);
	srand(time(NULL));
	_dxbase = (_xmax-_xmin)/(_numxNodes-1);
	_dybase = (_ymax - _ymin)/(_numyNodes-1);
	_nodeNumber = _numxNodes * _numyNodes;
	
	clock_t t1,t2;
	t1=clock();
	std::vector<Node> verticesbase(_nodeNumber);        // Intialize base vertices
	ofstream outdata;
	outdata.open (_filename);
	outdata << "Graph information: numx, numy, dx, dy\n";
	outdata << _numxNodes <<"\t"<<_numyNodes<<"\t"<<_xmin<<"\t"<<_xmax<<"\t"<<_ymin<<"\t"<<_ymax<<"\t"<<"\t"<<_radius<<"\n";
	outdata << "Number of nodes\n";
	outdata << _nodeNumber << "\n";
	outdata << "NodeId \t Coordinates \n";
	
	for( int a = 0; a < _nodeNumber; a++)
	{
		Node node =  gbase.addNode();
		verticesbase.at(a) = node;
		coordbase[node].x = _xmin + (a%(_numxNodes))*_dxbase;
		coordbase[node].y = _ymin + (floor(a/(_numxNodes)))*_dybase;
		outdata << a<< "\t" << coordbase[node].x << "\t" << coordbase[node].y << "\n";
	}	
	outdata << "Edge(i to  j)  Weight  \n";
	srand(time(NULL));
	
	double deltaxy = max(_dxbase,_dybase);
	int H = 0;
	while ( deltaxy*pow(2,H) <= _radius )
	{
		H++;
	}
	H = H - 1;
	
	for (int stage = 0 ; stage <= H; ++stage)
	{
		for (int i = 0; i < _numxNodes; i=i+pow(2,stage)){
		for (int j = 0; j < _numyNodes; j=j+pow(2,stage)){
			Node v = verticesbase.at(i + j*_numxNodes);
			_height[v] = stage;
		}}	
	}
	
// 	int factor[] = {1,4,8,16,32};  
// 	double distances[] = {pow(15*dxbase,2), pow(30*dxbase,2),pow(50*dxbase,2),pow(80*dxbase,2),pow(120*dxbase,2)};
// 	int factorSize = (sizeof(factor)/sizeof(*factor));
	
	double radius_square = _radius * _radius;
	int counter = 0;
	for (int a = 0; a< _nodeNumber; ++a){
		for (int b = a+1; b< _nodeNumber; ++b)
		{
			Node i = verticesbase.at(a), j = verticesbase.at(b);
			double x1 = coordbase[i].x, y1 = coordbase[i].y;
			double x2 = coordbase[j].x, y2 = coordbase[j].y;
			double dist = (x1-x2)*(x1-x2)+(y1-y2)*(y1-y2);
// 			for (int ik=0; ik< factorSize; ik++ ){ //TODO a%factor[ik]==0 && b%factor[ik]==0 && 
			if (dist < radius_square) // && ((double)rand()/(double)RAND_MAX) < 0.20 + 0.01*factor[ik]
			{
				Arc e = gbase.addArc(i, j); Arc em = gbase.addArc(j, i);                
				double arclength, arclengthv2;
				getGroundDistance(x1,y1,x2,y2,_numxNodes,_numyNodes,_dxbase, 
						_dybase, arclength,arclengthv2,_xmin,_ymin);
				lenbase.set(e, arclength);  lenbase.set(em, arclengthv2);		
				outdata << a << "\t" << b << "\t" << lenbase[e] <<"\n";
				outdata << b << "\t" << a << "\t" << lenbase[em] <<"\n";
                counter++;
			}
	}}
	
	t2=clock();	float diff ((float)t2-(float)t1);
	cout<<"|V(G)| = " << _numxNodes*_numyNodes << ", |E(G)| = "<< counter << ", running time: " <<diff/(CLOCKS_PER_SEC*60) << " minutes \n";
	outdata << "The shortest path is \n";
	outdata.close();
}


