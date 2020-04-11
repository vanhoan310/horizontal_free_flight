#include <math.h>
#include <iomanip>
#include <vector>
#include <algorithm>
#include <assert.h>
#include <lemon/list_graph.h>
#include <lemon/adaptors.h>  // filter for graph
#include <lemon/dijkstra.h>
#include <lemon/dim2.h>
#include <unordered_map>
#include <fstream>
#include <cstdlib>
#include <sstream>
#include "adaptivedijkstra33.h"
#include "adaptivedijkstra33.cpp"  

double getWindNorth(double x,double y) 
{ 
    return ((-2000-x)*exp(-(pow(x+2000,2) + pow(y-7000,2))/pow(10,8))-1.2*(8000-x)*exp(-(pow(x-8000,2) + pow(y+6000,2))/pow(10,7.5)))/24;
} 

double getWindEast(double x,double y) 
{ 
    return ((y-7000)*exp(-(pow(x+2000,2) + pow(y-7000,2))/pow(10,8))-1.2*(y+6000)*exp(-(pow(x-8000,2) + pow(y+6000,2))/pow(10,7.5)))/24;
}

const double distance2Points(const double& x1, const double& y1, const double& x2, const double& y2)
{
    return sqrt((x1-x2)*(x1-x2)+(y1-y2)*(y1-y2));
}

const double Line(const double& x1, const double& y1, const double& x2, const double& y2, const double& a)
{
    return (y2-y1)*(a-x1)/(x2-x1)+y1;
}

const int coordinate2VertexIndex(int numxbase,int numybase,double  dxbase,double dybase,double getx,double gety)
{
	const int n = round(getx/dxbase);
	const int m = round(gety/dybase);
	int result = (m*numxbase+n);
	return result;
}

inline int ipow(int base, int exp)
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

void getGroundDistancefromSource(const double& x1,const double& y1,const double& x2, 
const double& y2, const double& dx,const double& dy ,double & arclength)
{
	double minx = (x1<x2)? x1:x2, maxx = (x1<x2)? x2:x1;
	double miny = (y1<y2)? y1:y2, maxy = (y1<y2)? y2:y1; 
	double wCross, wTrack;  
	double vPlane = 850;
	double a,b,temp;
	double distance2P = distance2Points(x1,y1,x2,y2);
	arclength = 0.0; 
	if (distance2P == 0) return;
	double cosAngle = (y2-y1)/distance2P, sinAngle = (x2-x1)/distance2P;
		
	if (x1==x2 || y1==y2)
	{
		arclength = distance2P;
		return;
	}   
	vector<double> v1, v2, vectorMerged;
	int xIter = floor((maxx-minx)/dx);
	int yIter = floor((maxy-miny)/dy);
	if (x2 < x1)
	{
		for (int i=0;i<=xIter;++i) v1.push_back(minx+i*dx);
		v1.push_back(x1);
	}
	else
	{
		v1.push_back(x1);
		for (int i=xIter;i >=0;--i) v1.push_back(maxx - i*dx);
	}
	
	a = (x2-x1)/(y2-y1); b = x1 - a*y1;
	if (a>0)
	{
		if (y2 < y1)
		{
			for (int i=0;i<=yIter;++i) v2.push_back(a*(miny+i*dy)+b);
			v2.push_back(a*(y1)+b);
		}
		else
		{
			v2.push_back(a*(y1)+b);
			for (int i=yIter;i>=0;--i) v2.push_back(a*(maxy-i*dy)+b);
		}
	} 
	else
	{
		if (y2 < y1)
		{
			v2.push_back(a*(y1)+b);
			for (int i=yIter;i>=0;--i) v2.push_back(a*(miny+i*dy)+b);			
		}
		else
		{
			for (int i=0;i<=yIter;++i) v2.push_back(a*(maxy - i*dy)+b);
			v2.push_back(a*(y1)+b);			
		}
		
	}
	int vectorsize = v1.size() + v2.size();
	vectorMerged.resize(vectorsize);
	merge(v1.begin(), v1.end(), v2.begin(),v2.end(), vectorMerged.begin());
	double px1,px2,py1,py2;
	arclength = 0.0; 
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
			std::cout.precision(2);
// 			cout << arclength << "\t";
		}
// 		cout << px1 << "\t";
	}
	arclength = arclength*vPlane;
// 		cout << endl;
    return;
}

void getGroundDistancefromSink(const double& x1,const double& y1,const double& x2, 
const double& y2, const double& dx,const double& dy ,double & arclength)
{
	double minx = (x1<x2)? x1:x2, maxx = (x1<x2)? x2:x1;
	double miny = (y1<y2)? y1:y2, maxy = (y1<y2)? y2:y1; 
	double wCross, wTrack;  
	double vPlane = 850;
	double a,b,temp;
	double distance2P = distance2Points(x1,y1,x2,y2);
	arclength = 0.0; 
	if (distance2P == 0) return;
	double cosAngle = (y1-y2)/distance2P, sinAngle = (x1-x2)/distance2P;
		
	if (x1==x2 || y1==y2)
	{
		arclength = distance2P;
		return;
	}   
	vector<double> v1, v2, vectorMerged;
	int xIter = floor((maxx-minx)/dx);
	int yIter = floor((maxy-miny)/dy);
	if (x2 < x1)
	{
		for (int i=0;i<=xIter;++i) v1.push_back(minx+i*dx);
		v1.push_back(x1);
	}
	else
	{
		v1.push_back(x1);
		for (int i=xIter;i >=0;--i) v1.push_back(maxx - i*dx);
	}
	
	a = (x2-x1)/(y2-y1); b = x1 - a*y1;
	if (a>0)
	{
		if (y2 < y1)
		{
			for (int i=0;i<=yIter;++i) v2.push_back(a*(miny+i*dy)+b);
			v2.push_back(a*(y1)+b);
		}
		else
		{
			v2.push_back(a*(y1)+b);
			for (int i=yIter;i>=0;--i) v2.push_back(a*(maxy-i*dy)+b);
		}
	} 
	else
	{
		if (y2 < y1)
		{
			v2.push_back(a*(y1)+b);
			for (int i=yIter;i>=0;--i) v2.push_back(a*(miny+i*dy)+b);			
		}
		else
		{
			for (int i=0;i<=yIter;++i) v2.push_back(a*(maxy - i*dy)+b);
			v2.push_back(a*(y1)+b);			
		}
		
	}
	int vectorsize = v1.size() + v2.size();
	vectorMerged.resize(vectorsize);
	merge(v1.begin(), v1.end(), v2.begin(),v2.end(), vectorMerged.begin());
	double px1,px2,py1,py2;
	arclength = 0.0; 
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
			std::cout.precision(2);
// 			cout << arclength << "\t";
		}
// 		cout << px1 << "\t";
	}
	arclength = arclength*vPlane;
// 		cout << endl;
    return;
}
void runAlgorithm(string filename)
{	
	bool printStatistic = true;
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
	//srand(time(NULL));
	/* ----------- Compute the shortest path for base graph ------------------------*/
	int numxbase, numybase;
	double xmin, xmax, ymin,ymax;
	double dxbase, dybase;
	double radius;
// 	double fromPointx,fromPointy,toPointx,toPointy;
	int nodeNumber = 0;
	std::ifstream infile(filename); std::string line;
	if(printStatistic) 
        cout << "========================================================================================================\n";
	cout << "Start reading data ... \n";   
	std::getline(infile, line); // skip line 1
	std::getline(infile, line);
	std::istringstream iss(line); // read number of nodes
	
	iss >> numxbase>>numybase>>xmin>>xmax>>ymin>>ymax>>radius;
	dxbase =  (xmax-xmin)/(numxbase-1); 	dybase =(ymax - ymin)/(numybase-1);
	std::getline(infile, line); // skip line 3
	std::getline(infile, line);
	std::istringstream iss2(line); // read number of nodes
	iss2 >> nodeNumber;
	std::vector<Node> verticesbase(nodeNumber);        // Intialize base vertices
	std::getline(infile, line);
	
	// Read the coordinate from file
	int getnode;
	double getx, gety;
	while (std::getline(infile, line))
	{
		std::istringstream iss(line);
		if (!(iss >> getnode >> getx >> gety))
		{
			break; 
			
		} // error
		verticesbase.at(getnode) = gbase.addNode();
		coordbase[verticesbase.at(getnode)].x = getx;
		coordbase[verticesbase.at(getnode)].y = gety;
	}
	
	// Read the weight from file
	int edgeCount =0;
	int firstNode, secondNode; double weight;
	while (std::getline(infile, line))
	{
		std::istringstream iss(line);
		if (!(iss >> firstNode >>secondNode >> weight)) { break; } // error
		Node i = verticesbase.at(firstNode), j = verticesbase.at(secondNode);
		Arc e = gbase.addArc(i, j);
		lenbase.set(e, weight);
		edgeCount++;
	}
    cout << " ... finish reading data! \n";	
    // Input source and target 
	double fromPointxvec[] = {15000.0};//
	double fromPointyvec[] = {6000.0};//
	
	double toPointxvec[] =   {-7000.0};//
	double toPointyvec[] =   {5000.0};//
	
	int iterSize = (sizeof(fromPointxvec)/sizeof(*fromPointxvec));
	
    int stageCount = 0;
	for (int numberofStage = 4; numberofStage < 7; ++ numberofStage) // 3-7
	{	
		double astarFactor = 1;	
		std::vector<double> time1,time2;
		double timeCount = 0;
		double averageTime = 0, averageError = 0;
		for (int iter = 0; iter < iterSize; iter++ ) //TODO
		{
			// Add extra node from source
			verticesbase.push_back(gbase.addNode());
			coordbase[verticesbase.back()].x = fromPointxvec[iter];
			coordbase[verticesbase.back()].y = fromPointyvec[iter];			
			Node s = verticesbase.back();
			// Add extra node from sink
			verticesbase.push_back(gbase.addNode());
			coordbase[verticesbase.back()].x = toPointxvec[iter];
			coordbase[verticesbase.back()].y = toPointyvec[iter];			
			Node t = verticesbase.back();
			
			for (Graph::NodeIt v(gbase); v != INVALID; ++v)
			{
				// Add edge for source
				if(abs(coordbase[v].x - coordbase[s].x) < radius && 
					abs(coordbase[v].y - coordbase[s].y) < radius)
				{
					Arc e = gbase.addArc(s, v);
					getGroundDistancefromSource(coordbase[s].x, coordbase[s].y,coordbase[v].x, coordbase[v].y, 
					dxbase,dybase, weight);
					lenbase.set(e, weight);
					edgeCount++;
				}
			
				// Add edge for sink
				if(abs(coordbase[v].x - coordbase[t].x) < radius && 
					abs(coordbase[v].y - coordbase[t].y) < radius)
				{
					Arc e = gbase.addArc(v, t);
					getGroundDistancefromSink(coordbase[t].x, coordbase[t].y,coordbase[v].x, coordbase[v].y, 
					dxbase,dybase, weight);
					lenbase.set(e, weight);
					edgeCount++;
				}
				
			}
			Path<Graph> p; double d;
            clock_t t1,t2;
			if(printStatistic)cout << "Dijkstra algorithm\n";	
            
            t1=clock();	  
            dijkstra(gbase,lenbase).path(p).dist(d).run(s,t);
            t2=clock();    
            float diff ((float)t2-(float)t1);
            if(printStatistic)cout <<"shortest distance: "<< std::setprecision(5)<< std::fixed << d << ", time: "
                <<std::setprecision(3)<< std::fixed<< diff/(CLOCKS_PER_SEC)<< " (s).\n"; 
            timeCount = diff/(CLOCKS_PER_SEC);
            time1.push_back(timeCount);	//cout <<endl;
			
			// AdaptDijkstra
			vector<Node> pathNode;
			if(printStatistic) cout << "Adaptive Dijkstra algorithm\n";
			AdaptDijkstra test_alg1(gbase,lenbase,numberofStage, coordbase,numxbase,numybase); 
            
			test_alg1.run(s,t);
			time2.push_back(test_alg1.getRunTime());
			//Print shortest path
			ostringstream fileNameStream;                    // let this be empty
			fileNameStream << "runme" <<numberofStage-2 << iter << ".m"; // and pass "dice_" here
			string fileName = fileNameStream.str();  
			ofstream outdata;
			outdata.open(fileName.c_str());
			test_alg1.printShortestPath(outdata);
			
			if(printStatistic) 
			{	
                cout <<"shortest distance: "<< std::setprecision(5)<< std::fixed << test_alg1.getShortestDistance() 
                        << ", time: " << time2.at(iter) << " (s).\n";
                cout << "------------------------------------------------------\n";
				cout << "SUMMARY: |V| = "<< countNodes(gbase)<< ", |E| = "<<countArcs(gbase)<< ",s=(" <<std::setprecision(2)<< std::fixed<<coordbase[s].x<<
				", "<<coordbase[s].y<<")" << ",t=("<<coordbase[t].x<<
				", "<<coordbase[t].y<<")" << ", relative error: " << std::scientific << ((test_alg1.getShortestDistance()-d)/d)*100;
				cout  <<std::setprecision(2)<< " %, speed up: " <<std::setprecision(2)<<std::fixed<< timeCount/time2.at(iter)<< " times." << endl;
			}
			averageError += ((test_alg1.getShortestDistance()-d)/d)*100;
			averageTime += timeCount/time2.at(iter);
// 			cout << endl;
			
		}	
		
        cout << "K = " << numberofStage - 2 << ", astarFactor = " <<std::setprecision(2)<<std::fixed<< astarFactor << " ,average time speed up: " <<std::setprecision(3)<<std::fixed<< averageTime/iterSize << ", average relative error: "<< std::scientific << averageError/iterSize << " %" <<endl;
        cout <<endl;
        cout << "****************************************************************************************************\n";
        stageCount++;
	}
	
	infile.close();
	infile.clear();
}
int main(int argc, char* argv[])
{    
	// check arguments
    if (argc > 3 || argc < 2) {
        cout << "Wrong argument!\n";
        cout << "E.g: ./adaptiveDijkstra -gendata filename.txt: to generate a sample data\n";
        cout << "     ./adaptiveDijkstra filename.txt         : to run algorithm on filename.txt\n";
    } 
    else if (argc == 3){
        clock_t t1,t2;
        t1=clock();
        std::string filename = argv[2];
        // Generate a simple dataset 
        cout << "Start generating an instance ... \n";
    	int    numxNodes 	=  513; 	int numyNodes = 257;
    	double xmin = -20000.0, xmax = 20000.0, ymin = -10000.0, ymax = 10000.0;
    	double radius = 1280.0;
    	MakeGraph datafile(filename, numxNodes, numyNodes, xmin, xmax, ymin, ymax, radius);
    	datafile.makeGraphData();

        t2=clock();	float diff ((float)t2-(float)t1);	
        cout<<"Total runtime: " <<diff/(CLOCKS_PER_SEC*60) << " minutes\n";
    }
    else 
    {
        clock_t t1,t2;
        t1=clock();
        std::string filename = argv[1];
        // Run the algorithm 
        runAlgorithm(filename);
        
        t2=clock();	float diff ((float)t2-(float)t1);	
        cout<<"Total runtime: " <<diff/(CLOCKS_PER_SEC*60) << " minutes\n";
    }
}
