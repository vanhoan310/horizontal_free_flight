#ifndef TEST_ALGORITHM
#define TEST_ALGORITHM

using namespace lemon;
using namespace std;

class AdaptDijkstra {
	private:
	typedef ListDigraph Graph;
	typedef Graph::Node Node;
	typedef Graph::Arc Arc;
	typedef Graph::ArcMap<double> LengthMap;
	typedef Graph::ArcIt ArcIt;
	typedef Graph::NodeIt NodeIt;
	typedef lemon::dim2::Point<double> Point;
		
	Graph::NodeMap<Point> & _coordinate;	
	//Reference to the underlying digraph.
	const Graph &_graph;
	//Reference to the length map.
	const LengthMap &_length;
	// Distance map from source
	Graph::NodeMap<double> _distSource;
	// Distance map from sink
	Graph::NodeMap<double> _distSink;
	// Garbage variables
	Graph::NodeMap<double> potential;
	//Graph::NodeMap<bool>  processed;
	// Filter for subgraph
	Graph::NodeMap<bool> _filter;
	Graph::NodeMap<int> _height;
	// Store neighboorhood of vertices graph G for the
	Graph::NodeMap<vector<Node>> _neighbor;
	Graph::NodeMap<vector<Node>> _parent;
	Path<ListDigraph> _shortest_path;
	vector<Node> _path_Nodes;
	vector<Node> _old_path_Nodes;
	Node _source, _sink;
	double _astarFactor = 0;
	double _adaptiveFactor = 0;
	double _old_shortest_distance;	
	double _shortest_distance;
	int _stage, _maxStage;          
	const double _INFINITY = 100000000.0; 
	double _runtime;	
	double heuristicBallRadiusSquare = 1;
	
	const int _nodeNumber, _numxbase, _numybase;
	public:
	AdaptDijkstra (Graph & g, LengthMap& length, int numberofStage, Graph::NodeMap<Point> & coordinate,int numxbase, int numybase);
	// Compute power of integer: base^exp 
	inline int ipow(int base, int exp);
	// Compute heuristic function
	inline double calculateHeuristic (ListDigraph::NodeMap<double> &dist,  ListDigraph::NodeMap<vector<Node>> & _parent, ListDigraph::Node & v);
		
	// Find neighboorhood for the next stage
	void findNeighborNextStage(int currentStage);
	void findNeighborNextStageSmall(int currentStage);
	void constructParentNodes();
	
	void init (Node s, Node t); // Initialize the data
	
	void start (Node s, Node t);  // Start the alg in the first stage

	double getShortestDistance();
	//     double /*_old_shortest_distance*/;
	double _adaptFactor;
	
	vector<Node> getShortestPath();

	void processNextStage(Node s, Node t, double & time, bool savePathNodes, Graph::NodeMap<vector<Node>> & neighborSet);  // Process the next stage of AdaptDijkstra.
	void processNextStageBeforeLast(Node s, Node t, double & time, bool savePathNodes); 
	// Run the algorithm.
	double calculateDistance(ListDigraph::NodeMap<lemon::dim2::Point<double>>& coordbase,
	Node currentNode, Node sinkNode);
	void findLocalShortestPath(Node s, Node t);
	void run(Node s,Node t);
	void printShortestPath(std::ofstream& outdata);
	double getRunTime();
};

class MakeGraph {
private:
	typedef ListDigraph Graph;
	typedef Graph::Node Node;
	typedef Graph::Arc Arc;
	typedef Graph::ArcMap<double> LengthMap;
	typedef Graph::ArcIt ArcIt;
	typedef Graph::NodeIt NodeIt;	
	std::string _filename;
	int    _numxNodes 	=  513; 	int _numyNodes = 257;
	double _xmin = -20000.0, _xmax = 20000.0, _ymin = -10000.0, _ymax = 10000.0;
	double _dxbase 	=  0.0, _dybase = 0.0, _radius = 0.0;
	int    _nodeNumber = 0;
public:	
	MakeGraph (std::string filename, const int numxNodes, const int numyNodes, const double xmin, const double xmax, const double ymin, const double ymax, double radius);
	
	void makeGraphData();
	void getGroundDistance(const double& x1,const double& y1,const double& x2, 
		const double& y2,const int& numx,const int& numy,const double& dx,const double& dy
		,double & arclength, double & arclengthv2,double _xmin, double _ymin);
	//v
	double getWindNorth(double x,double y);
	//u
	double getWindEast(double x,double y);
	// Merge two increasing sorted vectors into one sorted vector in the increasing order
	std::vector<double> merge2Sorted( const std::vector<double>& left, const std::vector<double>& right );
	// Evaluate the line defined by two points (x1, y1) and (x2,y2) at a.
	const double Line(const double& x1, const double& y1, const double& x2, const double& y2, const double& a);
	// Compute distance between (x1,y1) and (x2,y2)
	const double distance2Points(const double& x1, const double& y1, const double& x2, const double& y2);
	
	
};
#endif
