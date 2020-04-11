// Adaptive dijkstra run with min du + c(u,v) + pv
#ifndef LEMON_ASTAR2DIR_H
#define LEMON_ASTAR2DIR_H

///\ingroup shortest_path
///\file
///\brief Astar2dir algorithm.

#include <limits>
#include <lemon/list_graph.h>
#include <lemon/bin_heap.h>
#include <lemon/bits/path_dump.h>
#include <lemon/core.h>
#include <lemon/error.h>
#include <lemon/maps.h>
#include <lemon/path.h>
#include <lemon/dim2.h>
#include <iostream>
using namespace std;
namespace lemon {

  /// \brief Default operation traits for the Astar2dir algorithm class.
  ///
  /// This operation traits class defines all computational operations and
  /// constants which are used in the Astar2dir algorithm.
  template <typename V>
  struct Astar2dirDefaultOperationTraits {
    /// \e
    typedef V Value;
    /// \brief Gives back the zero value of the type.
    static Value zero() {
      return static_cast<Value>(0);
    }
    /// \brief Gives back the sum of the given two elements.
    static Value plus(const Value& left, const Value& right) {
      return left + right;
    }
    /// \brief Gives back true only if the first value is less than the second.
    static bool less(const Value& left, const Value& right) {
      return left < right;
    }
    static bool lessorequal(const Value& left, const Value& right) {
      return left <= right;
    }
  };

  ///Default traits class of Astar2dir class.

  ///Default traits class of Astar2dir class.
  ///\tparam GR The type of the digraph.
  ///\tparam LEN The type of the length map.
  template<typename GR, typename LEN>
  struct Astar2dirDefaultTraits
  {
    ///The type of the digraph the algorithm runs on.
    typedef GR Digraph;

    ///The type of the map that stores the arc lengths.

    ///The type of the map that stores the arc lengths.
    ///It must conform to the \ref concepts::ReadMap "ReadMap" concept.
    typedef LEN LengthMap;
    ///The type of the arc lengths.
    typedef typename LEN::Value Value;

    /// Operation traits for %Astar2dir algorithm.

    /// This class defines the operations that are used in the algorithm.
    /// \see Astar2dirDefaultOperationTraits
    typedef Astar2dirDefaultOperationTraits<Value> OperationTraits;

    /// The cross reference type used by the heap.

    /// The cross reference type used by the heap.
    /// Usually it is \c Digraph::NodeMap<int>.
    typedef typename Digraph::template NodeMap<int> HeapCrossRef;
    ///Instantiates a \c HeapCrossRef.

    ///This function instantiates a \ref HeapCrossRef.
    /// \param g is the digraph, to which we would like to define the
    /// \ref HeapCrossRef.
    static HeapCrossRef *createHeapCrossRef(const Digraph &g)
    {
      return new HeapCrossRef(g);
    }

    typedef BinHeap<typename LEN::Value, HeapCrossRef, std::less<Value> > Heap;
    ///Instantiates a \c Heap.

    ///This function instantiates a \ref Heap.
    static Heap *createHeap(HeapCrossRef& r)
    {
      return new Heap(r);
    }

    ///\brief The type of the map that stores the predecessor
    ///arcs of the shortest paths.
    ///
    ///The type of the map that stores the predecessor
    ///arcs of the shortest paths.
    ///It must conform to the \ref concepts::WriteMap "WriteMap" concept.
    typedef typename Digraph::template NodeMap<typename Digraph::Arc> PredMap;
    ///Instantiates a \c PredMap.

    ///This function instantiates a \ref PredMap.
    ///\param g is the digraph, to which we would like to define the
    ///\ref PredMap.
    static PredMap *createPredMap(const Digraph &g)
    {
      return new PredMap(g);
    }

    ///The type of the map that indicates which nodes are processed.

    ///The type of the map that indicates which nodes are processed.
    ///It must conform to the \ref concepts::WriteMap "WriteMap" concept.
    ///By default, it is a NullMap.
    typedef NullMap<typename Digraph::Node,bool> ProcessedMap;
    ///Instantiates a \c ProcessedMap.

    ///This function instantiates a \ref ProcessedMap.
    ///\param g is the digraph, to which
    ///we would like to define the \ref ProcessedMap.
#ifdef DOXYGEN
    static ProcessedMap *createProcessedMap(const Digraph &g)
#else
    static ProcessedMap *createProcessedMap(const Digraph &)
#endif
    {
      return new ProcessedMap();
    }

    ///The type of the map that stores the distances of the nodes.

    ///The type of the map that stores the distances of the nodes.
    ///It must conform to the \ref concepts::WriteMap "WriteMap" concept.
    typedef typename Digraph::template NodeMap<typename LEN::Value> DistMap;
    ///Instantiates a \c DistMap.

    ///This function instantiates a \ref DistMap.
    ///\param g is the digraph, to which we would like to define
    ///the \ref DistMap.
    static DistMap *createDistMap(const Digraph &g)
    {
      return new DistMap(g);
    }
  };

  ///%Astar2dir algorithm class.

  /// \ingroup shortest_path
  ///This class provides an efficient implementation of the %Astar2dir algorithm.
  ///
  ///The %Astar2dir algorithm solves the single-source shortest path problem
  ///when all arc lengths are non-negative. If there are negative lengths,
  ///the BellmanFord algorithm should be used instead.
  ///
  ///The arc lengths are passed to the algorithm using a
  ///\ref concepts::ReadMap "ReadMap",
  ///so it is easy to change it to any kind of length.
  ///The type of the length is determined by the
  ///\ref concepts::ReadMap::Value "Value" of the length map.
  ///It is also possible to change the underlying priority heap.
  ///
  ///There is also a \ref astar2dir() "function-type interface" for the
  ///%Astar2dir algorithm, which is convenient in the simplier cases and
  ///it can be used easier.
  ///
  ///\tparam GR The type of the digraph the algorithm runs on.
  ///The default type is \ref ListDigraph.
  ///\tparam LEN A \ref concepts::ReadMap "readable" arc map that specifies
  ///the lengths of the arcs.
  ///It is read once for each arc, so the map may involve in
  ///relatively time consuming process to compute the arc lengths if
  ///it is necessary. The default map type is \ref
  ///concepts::Digraph::ArcMap "GR::ArcMap<int>".
  ///\tparam TR The traits class that defines various types used by the
  ///algorithm. By default, it is \ref Astar2dirDefaultTraits
  ///"Astar2dirDefaultTraits<GR, LEN>".
  ///In most cases, this parameter should not be set directly,
  ///consider to use the named template parameters instead.
#ifdef DOXYGEN
  template <typename GR, typename LEN, typename TR>
#else
  template <typename GR=ListDigraph,
            typename LEN=typename GR::template ArcMap<int>,
            typename TR=Astar2dirDefaultTraits<GR,LEN> >
#endif
  class Astar2dir {
  public:

    ///The type of the digraph the algorithm runs on.
    typedef typename TR::Digraph Digraph;

    ///The type of the arc lengths.
    typedef typename TR::Value Value;
    ///The type of the map that stores the arc lengths.
    typedef typename TR::LengthMap LengthMap;
    ///\brief The type of the map that stores the predecessor arcs of the
    ///shortest paths.
    typedef typename TR::PredMap PredMap;
    ///The type of the map that stores the distances of the nodes.
    typedef typename TR::DistMap DistMap;
    ///The type of the map that indicates which nodes are processed.
    typedef typename TR::ProcessedMap ProcessedMap;
    ///The type of the paths.
    typedef PredMapPath<Digraph, PredMap> Path;
    ///The cross reference type used for the current heap.
    typedef typename TR::HeapCrossRef HeapCrossRef;
    ///The heap type used by the algorithm.
    typedef typename TR::Heap Heap;
    /// \brief The \ref lemon::Astar2dirDefaultOperationTraits
    /// "operation traits class" of the algorithm.
    typedef typename TR::OperationTraits OperationTraits;

    ///The \ref lemon::Astar2dirDefaultTraits "traits class" of the algorithm.
    typedef TR Traits;

  private:

    typedef typename Digraph::Node Node;
    typedef typename Digraph::NodeIt NodeIt;
    typedef typename Digraph::Arc Arc;
    typedef typename Digraph::OutArcIt OutArcIt;
    typedef typename Digraph::InArcIt InArcIt;

    //Pointer to the underlying digraph.
    const Digraph *G;
    //Pointer to the length map.
    const LengthMap *_length;
    //Pointer to the map of predecessors arcs.
    PredMap *_pred;
    //Indicates if _pred is locally allocated (true) or not.
    bool local_pred;
    //Pointer to the map of distances.
    DistMap *_dist;
    //Indicates if _dist is locally allocated (true) or not.
    bool local_dist;
    //Pointer to the map of processed status of the nodes.
    ProcessedMap *_processed;
    //Indicates if _processed is locally allocated (true) or not.
    bool local_processed;
    //Pointer to the heap cross references.
    HeapCrossRef *_heap_cross_ref;
    //Indicates if _heap_cross_ref is locally allocated (true) or not.
    bool local_heap_cross_ref;
    //Pointer to the heap.
    Heap *_heap;
    //Indicates if _heap is locally allocated (true) or not.
    bool local_heap;
    int _exploreEdges;
    Value _shortest_distance;
    //TODO
    //ListDigraph::NodeMap<lemon::dim2::Point<double>> _coordbase;
    
    // Heuristic function for wind effect. 
   
    //Creates the maps if necessary.
    void create_maps()
    {
      if(!_pred) {
        local_pred = true;
        _pred = Traits::createPredMap(*G);
      }
      if(!_dist) {
        local_dist = true;
        _dist = Traits::createDistMap(*G);
      }
      if(!_processed) {
        local_processed = true;
        _processed = Traits::createProcessedMap(*G);
      }
      if (!_heap_cross_ref) {
        local_heap_cross_ref = true;
        _heap_cross_ref = Traits::createHeapCrossRef(*G);
      }
      if (!_heap) {
        local_heap = true;
        _heap = Traits::createHeap(*_heap_cross_ref);
      }
    }

  public:

    typedef Astar2dir Create;

    ///\name Named Template Parameters

    ///@{

    template <class T>
    struct SetPredMapTraits : public Traits {
      typedef T PredMap;
      static PredMap *createPredMap(const Digraph &)
      {
        LEMON_ASSERT(false, "PredMap is not initialized");
        return 0; // ignore warnings
      }
    };
    ///\brief \ref named-templ-param "Named parameter" for setting
    ///\c PredMap type.
    ///
    ///\ref named-templ-param "Named parameter" for setting
    ///\c PredMap type.
    ///It must conform to the \ref concepts::WriteMap "WriteMap" concept.
    template <class T>
    struct SetPredMap
      : public Astar2dir< Digraph, LengthMap, SetPredMapTraits<T> > {
      typedef Astar2dir< Digraph, LengthMap, SetPredMapTraits<T> > Create;
    };

    template <class T>
    struct SetDistMapTraits : public Traits {
      typedef T DistMap;
      static DistMap *createDistMap(const Digraph &)
      {
        LEMON_ASSERT(false, "DistMap is not initialized");
        return 0; // ignore warnings
      }
    };
    ///\brief \ref named-templ-param "Named parameter" for setting
    ///\c DistMap type.
    ///
    ///\ref named-templ-param "Named parameter" for setting
    ///\c DistMap type.
    ///It must conform to the \ref concepts::WriteMap "WriteMap" concept.
    template <class T>
    struct SetDistMap
      : public Astar2dir< Digraph, LengthMap, SetDistMapTraits<T> > {
      typedef Astar2dir< Digraph, LengthMap, SetDistMapTraits<T> > Create;
    };

    template <class T>
    struct SetProcessedMapTraits : public Traits {
      typedef T ProcessedMap;
      static ProcessedMap *createProcessedMap(const Digraph &)
      {
        LEMON_ASSERT(false, "ProcessedMap is not initialized");
        return 0; // ignore warnings
      }
    };
    ///\brief \ref named-templ-param "Named parameter" for setting
    ///\c ProcessedMap type.
    ///
    ///\ref named-templ-param "Named parameter" for setting
    ///\c ProcessedMap type.
    ///It must conform to the \ref concepts::WriteMap "WriteMap" concept.
    template <class T>
    struct SetProcessedMap
      : public Astar2dir< Digraph, LengthMap, SetProcessedMapTraits<T> > {
      typedef Astar2dir< Digraph, LengthMap, SetProcessedMapTraits<T> > Create;
    };

    struct SetStandardProcessedMapTraits : public Traits {
      typedef typename Digraph::template NodeMap<bool> ProcessedMap;
      static ProcessedMap *createProcessedMap(const Digraph &g)
      {
        return new ProcessedMap(g);
      }
    };
    ///\brief \ref named-templ-param "Named parameter" for setting
    ///\c ProcessedMap type to be <tt>Digraph::NodeMap<bool></tt>.
    ///
    ///\ref named-templ-param "Named parameter" for setting
    ///\c ProcessedMap type to be <tt>Digraph::NodeMap<bool></tt>.
    ///If you don't set it explicitly, it will be automatically allocated.
    struct SetStandardProcessedMap
      : public Astar2dir< Digraph, LengthMap, SetStandardProcessedMapTraits > {
      typedef Astar2dir< Digraph, LengthMap, SetStandardProcessedMapTraits >
      Create;
    };

    template <class H, class CR>
    struct SetHeapTraits : public Traits {
      typedef CR HeapCrossRef;
      typedef H Heap;
      static HeapCrossRef *createHeapCrossRef(const Digraph &) {
        LEMON_ASSERT(false, "HeapCrossRef is not initialized");
        return 0; // ignore warnings
      }
      static Heap *createHeap(HeapCrossRef &)
      {
        LEMON_ASSERT(false, "Heap is not initialized");
        return 0; // ignore warnings
      }
    };
    ///\brief \ref named-templ-param "Named parameter" for setting
    ///heap and cross reference types
    ///
    ///\ref named-templ-param "Named parameter" for setting heap and cross
    ///reference types. If this named parameter is used, then external
    ///heap and cross reference objects must be passed to the algorithm
    ///using the \ref heap() function before calling \ref run(Node) "run()"
    ///or \ref init().
    ///\sa SetStandardHeap
    template <class H, class CR = typename Digraph::template NodeMap<int> >
    struct SetHeap
      : public Astar2dir< Digraph, LengthMap, SetHeapTraits<H, CR> > {
      typedef Astar2dir< Digraph, LengthMap, SetHeapTraits<H, CR> > Create;
    };

    template <class H, class CR>
    struct SetStandardHeapTraits : public Traits {
      typedef CR HeapCrossRef;
      typedef H Heap;
      static HeapCrossRef *createHeapCrossRef(const Digraph &G) {
        return new HeapCrossRef(G);
      }
      static Heap *createHeap(HeapCrossRef &R)
      {
        return new Heap(R);
      }
    };
    ///\brief \ref named-templ-param "Named parameter" for setting
    ///heap and cross reference types with automatic allocation
    ///
    ///\ref named-templ-param "Named parameter" for setting heap and cross
    ///reference types with automatic allocation.
    ///They should have standard constructor interfaces to be able to
    ///automatically created by the algorithm (i.e. the digraph should be
    ///passed to the constructor of the cross reference and the cross
    ///reference should be passed to the constructor of the heap).
    ///However, external heap and cross reference objects could also be
    ///passed to the algorithm using the \ref heap() function before
    ///calling \ref run(Node) "run()" or \ref init().
    ///\sa SetHeap
    template <class H, class CR = typename Digraph::template NodeMap<int> >
    struct SetStandardHeap
      : public Astar2dir< Digraph, LengthMap, SetStandardHeapTraits<H, CR> > {
      typedef Astar2dir< Digraph, LengthMap, SetStandardHeapTraits<H, CR> >
      Create;
    };

    template <class T>
    struct SetOperationTraitsTraits : public Traits {
      typedef T OperationTraits;
    };

    /// \brief \ref named-templ-param "Named parameter" for setting
    ///\c OperationTraits type
    ///
    ///\ref named-templ-param "Named parameter" for setting
    ///\c OperationTraits type.
    /// For more information, see \ref Astar2dirDefaultOperationTraits.
    template <class T>
    struct SetOperationTraits
      : public Astar2dir<Digraph, LengthMap, SetOperationTraitsTraits<T> > {
      typedef Astar2dir<Digraph, LengthMap, SetOperationTraitsTraits<T> >
      Create;
    };

    ///@}

  protected:

    Astar2dir() {}

  public:

    ///Constructor.

    ///Constructor.
    ///\param g The digraph the algorithm runs on.
    ///\param length The length map used by the algorithm.
    Astar2dir(const Digraph& g, const LengthMap& length, bool reverse=false) :
      G(&g), _length(&length),
      _pred(NULL), local_pred(false),
      _dist(NULL), local_dist(false),
      _processed(NULL), local_processed(false),
      _heap_cross_ref(NULL), local_heap_cross_ref(false),
      _heap(NULL), local_heap(false), _exploreEdges(0)
    { }

    ///Destructor.
    ~Astar2dir()
    {
      if(local_pred) delete _pred;
      if(local_dist) delete _dist;
      if(local_processed) delete _processed;
      if(local_heap_cross_ref) delete _heap_cross_ref;
      if(local_heap) delete _heap;
    }

    ///Sets the length map.

    ///Sets the length map.
    ///\return <tt> (*this) </tt>
    
    Astar2dir &lengthMap(const LengthMap &m)
    {
      _length = &m;
      return *this;
    }
    
    ///Sets the map that stores the predecessor arcs.

    ///Sets the map that stores the predecessor arcs.
    ///If you don't use this function before calling \ref run(Node) "run()"
    ///or \ref init(), an instance will be allocated automatically.
    ///The destructor deallocates this automatically allocated map,
    ///of course.
    ///\return <tt> (*this) </tt>
    Astar2dir &predMap(PredMap &m)
    {
      if(local_pred) {
        delete _pred;
        local_pred=false;
      }
      _pred = &m;
      return *this;
    }

    ///Sets the map that indicates which nodes are processed.

    ///Sets the map that indicates which nodes are processed.
    ///If you don't use this function before calling \ref run(Node) "run()"
    ///or \ref init(), an instance will be allocated automatically.
    ///The destructor deallocates this automatically allocated map,
    ///of course.
    ///\return <tt> (*this) </tt>
    Astar2dir &processedMap(ProcessedMap &m)
    {
      if(local_processed) {
        delete _processed;
        local_processed=false;
      }
      _processed = &m;
      return *this;
    }

    ///Sets the map that stores the distances of the nodes.

    ///Sets the map that stores the distances of the nodes calculated by the
    ///algorithm.
    ///If you don't use this function before calling \ref run(Node) "run()"
    ///or \ref init(), an instance will be allocated automatically.
    ///The destructor deallocates this automatically allocated map,
    ///of course.
    ///\return <tt> (*this) </tt>
    Astar2dir &distMap(DistMap &m)
    {
      if(local_dist) {
        delete _dist;
        local_dist=false;
      }
      _dist = &m;
      return *this;
    }

    ///Sets the heap and the cross reference used by algorithm.

    ///Sets the heap and the cross reference used by algorithm.
    ///If you don't use this function before calling \ref run(Node) "run()"
    ///or \ref init(), heap and cross reference instances will be
    ///allocated automatically.
    ///The destructor deallocates these automatically allocated objects,
    ///of course.
    ///\return <tt> (*this) </tt>
    Astar2dir &heap(Heap& hp, HeapCrossRef &cr)
    {
      if(local_heap_cross_ref) {
        delete _heap_cross_ref;
        local_heap_cross_ref=false;
      }
      _heap_cross_ref = &cr;
      if(local_heap) {
        delete _heap;
        local_heap=false;
      }
      _heap = &hp;
      return *this;
    }

  private:

    void finalizeNodeData(Node v,Value dst)
    {
      _processed->set(v,true);
      //_dist->set(v, dst);
    }

  public:

    ///\name Execution Control
    ///The simplest way to execute the %Astar2dir algorithm is to use
    ///one of the member functions called \ref run(Node) "run()".\n
    ///If you need better control on the execution, you have to call
    ///\ref init() first, then you can add several source nodes with
    ///\ref addSource(). Finally the actual path computation can be
    ///performed with one of the \ref start() functions.

    ///@{

    ///\brief Initializes the internal data structures.
    ///
    ///Initializes the internal data structures.
    void init()
    {
      create_maps();
      _heap->clear();
      for ( NodeIt u(*G) ; u!=INVALID ; ++u ) {
        _pred->set(u,INVALID);
        _processed->set(u,false);
        _heap_cross_ref->set(u,Heap::PRE_HEAP);
      }
    }

    ///Adds a new source node.

    ///Adds a new source node to the priority heap.
    ///The optional second parameter is the initial distance of the node.
    ///
    ///The function checks if the node has already been added to the heap and
    ///it is pushed to the heap only if either it was not in the heap
    ///or the shortest path found till then is shorter than \c dst.
    void addSource(Node s,Value dst=OperationTraits::zero())
    {
      if(_heap->state(s) != Heap::IN_HEAP) {
        _heap->push(s,dst);
      } else if(OperationTraits::less((*_heap)[s], dst)) {
        _heap->set(s,dst);
        _pred->set(s,INVALID);
      }
    }

    ///Processes the next node in the priority heap

    ///Processes the next node in the priority heap.
    ///
    ///\return The processed node.
    ///
    ///\warning The priority heap must not be empty.
/*inline double calculateHeuristic(ListDigraph::NodeMap<lemon::dim2::Point<double>>& coordbase,
      Node currentNode, Node sinkNode, double factor
    ) 
    {
      const double x1 = coordbase[currentNode].x, y1 = coordbase[currentNode].y;
      const double x2 = coordbase[sinkNode].x,    y2 = coordbase[sinkNode].y;
      return sqrt((x1-x2)*(x1-x2)+(y1-y2)*(y1-y2))*factor;
    }*/
 
inline double calculateHeuristic (ListDigraph::NodeMap<double> &dist, ListDigraph::NodeMap<vector<Node>> & _neighbor, ListDigraph::Node & v, double _INFINITY, ListDigraph::NodeMap<vector<Node>> & _parent)
{
// 	double heuristicValue = dist[v];
// 	for (int i = 0; i< _neighbor[v].size(); ++i) 
// 	{  
// 	if (dist[_neighbor[v].at(i)] < heuristicValue) 
// 	heuristicValue = dist[_neighbor[v].at(i)]; 
// 	}
// 	return heuristicValue; //TODO 
//  	if(_parent[v].at(0) == lemon::INVALID) return dist[v]; 
	
	
	
// 	double heuristicValue = dist[v];
// 	for (int i = 0; i< _neighbor[v].size(); ++i) 
// 	{  
	if(_neighbor[v].empty()) return dist[v]; 
	double heuristicValue = dist[_neighbor[v].at(0)];
	for (int i = 1; i< _neighbor[v].size(); ++i) 
	{  
		Node w = _neighbor[v].at(i);
// 		if (dist[w] == _INFINITY) dist[w] = calculateHeuristic(dist,_parent,w,_INFINITY,_parent);
		if (dist[w] < heuristicValue )
			heuristicValue = dist[w];
	}
//  	dist[v] = heuristicValue; //TODO co nen gan luon gia tri cho thang nay ko ???
// 	if (heuristicValue == _INFINITY)
// 	{
// 		for (int i = 0; i< _parent[v].size(); ++i) 
// 		{  
// 			Node w = _parent[v].at(i);
// 			if (dist[w] < heuristicValue )
// 				heuristicValue = dist[w];
// 		}
// 	}
// 	if (heuristicValue == _INFINITY) return 2000.0;
	return heuristicValue; 
}  

bool runastar(Node& s, Node& t,ListDigraph::NodeMap<vector<Node>> & _neighbor,
		ListDigraph::NodeMap<double> & _dists,ListDigraph::NodeMap<double> & _distt,
		ListDigraph::NodeMap<double> & potential, double shortest_distance, double _INFINITY, ListDigraph::NodeMap<vector<Node>> & _parent, double adaptFactor)
{
	
	//for (NodeIt v(*G); v != INVALID; ++v){
	//	potential[v] = calculateHeuristic(_distt,_neighbor,v, _INFINITY);   }
// 	int count = 0;
 	for (NodeIt v(*G); v != INVALID; ++v){
 		potential[v] = -1.0; }
//  	for(OutArcIt e(*G,s); e!=INVALID; ++e) 
// 	{
// 		Node w=G->target(e);
// 		potential[w] = 0;
// 	}
// 	for(InArcIt e(*G,t); e!=INVALID; ++e) 
// 	{
// 		Node w=G->source(e);
// 		potential[w] = 0;
// 	}
	potential[s] = 0;  
	potential[t] = 0;
	
// 	shortest_distance *= 1.0001; // TODO EPS machine error
// 	shortest_distance *= adaptFactor; //1.06  1.02 good
	shortest_distance += 0.00001; //1.06  1.02 good
	init();
	addSource(s);  
	Value testvalue;
	while ( !_heap->empty() && _heap->top()!=t ) 
	{
		Node v=_heap->top();
		Value oldvalue=_heap->prio();
		_heap->pop();
		//finalizeNodeData(v,oldvalue);
		_dists[v] = oldvalue;
		for(OutArcIt e(*G,v); e!=INVALID; ++e) 
		{
			Node w=G->target(e);	
			if (potential[w] <= -1.0) 
			{
				potential[w] = calculateHeuristic(_distt,_neighbor,w,_INFINITY,_parent); 
			}
// 			cout << potential[w] << "\t";
			Value newvalue = OperationTraits::plus(oldvalue, (*_length)[e]);
			if (OperationTraits::lessorequal(OperationTraits::plus(newvalue, potential[w]), shortest_distance) )
			{
				switch(_heap->state(w)) {
				case Heap::PRE_HEAP:			
					_heap->push(w,newvalue);
					_pred->set(w,e);					
					break;
				case Heap::IN_HEAP:
				{
					if ( OperationTraits::less(newvalue, (*_heap)[w]) ) 
					{
						_heap->decrease(w, newvalue);
						_pred->set(w,e);
					}
				}
					break;
				case Heap::POST_HEAP:
					break;
				}
			}
		}
	}
		
	if ( !_heap->empty() ) 
	{
		Node v=_heap->top();
		_shortest_distance =_heap->prio();	  
		_heap->pop();
		_dists[v] = _shortest_distance;
	}

// 	for (NodeIt v(*G); v != INVALID; ++v){
// 		if (potential[v] > 0) _distt[v] = min(_dists[v], potential[v]);}
// 	potential[s] = 0;  
// 	potential[t] = 0;
	
	if((*_heap_cross_ref)[t] != Heap::POST_HEAP) 
	{
		cout << "There is no path from source to sink \n";
		_shortest_distance = 2*_INFINITY;
	}
	return (*_heap_cross_ref)[t] == Heap::POST_HEAP;
}


bool runastarAdaptive(Node& s, Node& t, ListDigraph::NodeMap<double> & _dists, double heuristicPercent, double _INFINITY)
{
	init();
	addSource(s);   
	int counter = 0;
	while ( !_heap->empty() && _heap->top()!=t ) 
	{
		Node v=_heap->top();
		Value oldvalue=_heap->prio();
		_heap->pop();
		//finalizeNodeData(v,oldvalue);
		_dists[v] = oldvalue;
		counter ++ ;
		for(OutArcIt e(*G,v); e!=INVALID; ++e) 
		{
			Node w=G->target(e);	    		
			switch(_heap->state(w)) 
			{
				case Heap::PRE_HEAP:
					_heap->push(w,OperationTraits::plus(oldvalue, (*_length)[e]));
					_pred->set(w,e);
					break;
				case Heap::IN_HEAP:
				{
					Value newvalue = OperationTraits::plus(oldvalue, (*_length)[e]);
					if ( OperationTraits::less(newvalue, (*_heap)[w]) ) {
					_heap->decrease(w, newvalue);
					_pred->set(w,e);
					}
				}
				break;
				case Heap::POST_HEAP:			
				break;
			}
		}
	}
	if ( !_heap->empty() ) 
	{
		Node v=_heap->top();
		_shortest_distance =_heap->prio();
		_dists[v] = _shortest_distance;
		_heap->pop();
	}
	if((*_heap_cross_ref)[t] != Heap::POST_HEAP) 
	{
		cout << "There is no path from source to sink in first stage \n";
		_shortest_distance = 2*_INFINITY;
	}
	else
	{
		int maxCounter = round(counter*heuristicPercent);
		//cout << maxCounter << "maxCounter";
		counter = 0;
		while ( !_heap->empty() && counter < maxCounter ) 
		{
			Node v=_heap->top();
			Value oldvalue=_heap->prio();
			_heap->pop();
			//finalizeNodeData(v,oldvalue);
			_dists[v] = oldvalue;
			counter ++ ;
			for(OutArcIt e(*G,v); e!=INVALID; ++e) 
			{
				Node w=G->target(e);	    		
				switch(_heap->state(w)) 
				{
					case Heap::PRE_HEAP:
						_heap->push(w,OperationTraits::plus(oldvalue, (*_length)[e]));
						_pred->set(w,e);
						break;
					case Heap::IN_HEAP:
					{
						Value newvalue = OperationTraits::plus(oldvalue, (*_length)[e]);
						if ( OperationTraits::less(newvalue, (*_heap)[w]) ) {
						_heap->decrease(w, newvalue);
						_pred->set(w,e);
						}
					}
					break;
					case Heap::POST_HEAP:			
					break;
					}
			}
		}
	}
	
	return (*_heap_cross_ref)[t] == Heap::POST_HEAP;
}

	
int getExploreEdges() const
{
	return _exploreEdges;
}  
	
	
    
    
    ///Returns the next node to be processed or \c INVALID if the
    ///priority heap is empty.
    Node nextNode() const
    {
      return !_heap->empty()?_heap->top():INVALID;
    }

    ///Returns \c false if there are nodes to be processed.

    ///Returns \c false if there are nodes to be processed
    ///in the priority heap.
    bool emptyQueue() const { return _heap->empty(); }

    ///Returns the number of the nodes to be processed.

    ///Returns the number of the nodes to be processed
    ///in the priority heap.
    int queueSize() const { return _heap->size(); }

    ///Executes the algorithm.

    ///Executes the algorithm.
    ///
    ///This method runs the %Astar2dir algorithm from the root node(s)
    ///in order to compute the shortest path to each node.
    ///
    ///The algorithm computes
    ///- the shortest path tree (forest),
    ///- the distance of each node from the root(s).
    ///
    ///\pre init() must be called and at least one root node should be
    ///added with addSource() before using this function.
    ///
    ///\note <tt>d.start()</tt> is just a shortcut of the following code.
    ///\code
    ///  while ( !d.emptyQueue() ) {
    ///    d.processNextNode();
    ///  }
    ///\endcode
    
    ///Executes the algorithm until the given target node is processed.

    ///Executes the algorithm until the given target node is processed.
    ///
    ///This method runs the %Astar2dir algorithm from the root node(s)
    ///in order to compute the shortest path to \c t.
    ///
    ///The algorithm computes
    ///- the shortest path to \c t,
    ///- the distance of \c t from the root(s).
    ///
    ///\pre init() must be called and at least one root node should be
    ///added with addSource() before using this function.
    
    ///Executes the algorithm until a condition is met.

    ///Executes the algorithm until a condition is met.
    ///
    ///This method runs the %Astar2dir algorithm from the root node(s) in
    ///order to compute the shortest path to a node \c v with
    /// <tt>nm[v]</tt> true, if such a node can be found.
    ///
    ///\param nm A \c bool (or convertible) node map. The algorithm
    ///will stop when it reaches a node \c v with <tt>nm[v]</tt> true.
    ///
    ///\return The reached node \c v with <tt>nm[v]</tt> true or
    ///\c INVALID if no such node was found.
    ///
    ///\pre init() must be called and at least one root node should be
    ///added with addSource() before using this function.
    
    ///Runs the algorithm from the given source node.

    ///This method runs the %Astar2dir algorithm from node \c s
    ///in order to compute the shortest path to each node.
    ///
    ///The algorithm computes
    ///- the shortest path tree,
    ///- the distance of each node from the root.
    ///
    ///\note <tt>d.run(s)</tt> is just a shortcut of the following code.
    ///\code
    ///  d.init();
    ///  d.addSource(s);
    ///  d.start();
    ///\endcode
   
    ///Finds the shortest path between \c s and \c t.

    ///This method runs the %Astar2dir algorithm from node \c s
    ///in order to compute the shortest path to node \c t
    ///(it stops searching when \c t is processed).
    ///
    ///\return \c true if \c t is reachable form \c s.
    ///
    ///\note Apart from the return value, <tt>d.run(s,t)</tt> is just a
    ///shortcut of the following code.
    ///\code
    ///  d.init();
    ///  d.addSource(s);
    ///  d.start(t);
    ///\endcode
    bool run(Node s,Node t) {
      init();
      addSource(s);
      start(t);
      return (*_heap_cross_ref)[t] == Heap::POST_HEAP;
    }

    ///@}

    ///\name Query Functions
    ///The results of the %Astar2dir algorithm can be obtained using these
    ///functions.\n
    ///Either \ref run(Node) "run()" or \ref init() should be called
    ///before using them.

    ///@{

    ///The shortest path to the given node.

    ///Returns the shortest path to the given node from the root(s).
    ///
    ///\warning \c t should be reached from the root(s).
    ///
    ///\pre Either \ref run(Node) "run()" or \ref init()
    ///must be called before using this function.
    Path path(Node t) const { return Path(*G, *_pred, t); }

    ///The distance of the given node from the root(s).

    ///Returns the distance of the given node from the root(s).
    ///
    ///\warning If node \c v is not reached from the root(s), then
    ///the return value of this function is undefined.
    ///
    ///\pre Either \ref run(Node) "run()" or \ref init()
    ///must be called before using this function.
    Value dist(Node v) const { return _shortest_distance; }

    ///\brief Returns the 'previous arc' of the shortest path tree for
    ///the given node.
    ///
    ///This function returns the 'previous arc' of the shortest path
    ///tree for the node \c v, i.e. it returns the last arc of a
    ///shortest path from a root to \c v. It is \c INVALID if \c v
    ///is not reached from the root(s) or if \c v is a root.
    ///
    ///The shortest path tree used here is equal to the shortest path
    ///tree used in \ref predNode() and \ref predMap().
    ///
    ///\pre Either \ref run(Node) "run()" or \ref init()
    ///must be called before using this function.
    Arc predArc(Node v) const { return (*_pred)[v]; }

    ///\brief Returns the 'previous node' of the shortest path tree for
    ///the given node.
    ///
    ///This function returns the 'previous node' of the shortest path
    ///tree for the node \c v, i.e. it returns the last but one node
    ///of a shortest path from a root to \c v. It is \c INVALID
    ///if \c v is not reached from the root(s) or if \c v is a root.
    ///
    ///The shortest path tree used here is equal to the shortest path
    ///tree used in \ref predArc() and \ref predMap().
    ///
    ///\pre Either \ref run(Node) "run()" or \ref init()
    ///must be called before using this function.
    Node predNode(Node v) const {  return (*_pred)[v]==INVALID ? INVALID:
                                  G->source((*_pred)[v]);
      }

    ///\brief Returns a const reference to the node map that stores the
    ///distances of the nodes.
    ///
    ///Returns a const reference to the node map that stores the distances
    ///of the nodes calculated by the algorithm.
    ///
    ///\pre Either \ref run(Node) "run()" or \ref init()
    ///must be called before using this function.
    const DistMap &distMap() const { return *_dist;}

    ///\brief Returns a const reference to the node map that stores the
    ///predecessor arcs.
    ///
    ///Returns a const reference to the node map that stores the predecessor
    ///arcs, which form the shortest path tree (forest).
    ///
    ///\pre Either \ref run(Node) "run()" or \ref init()
    ///must be called before using this function.
    const PredMap &predMap() const { return *_pred;}

    ///Checks if the given node is reached from the root(s).

    ///Returns \c true if \c v is reached from the root(s).
    ///
    ///\pre Either \ref run(Node) "run()" or \ref init()
    ///must be called before using this function.
    bool reached(Node v) const { return (*_heap_cross_ref)[v] !=
                                        Heap::PRE_HEAP; }

    ///Checks if a node is processed.

    ///Returns \c true if \c v is processed, i.e. the shortest
    ///path to \c v has already found.
    ///
    ///\pre Either \ref run(Node) "run()" or \ref init()
    ///must be called before using this function.
    bool processed(Node v) const { return (*_heap_cross_ref)[v] ==
                                          Heap::POST_HEAP; }

    ///The current distance of the given node from the root(s).

    ///Returns the current distance of the given node from the root(s).
    ///It may be decreased in the following processes.
    ///
    ///\pre Either \ref run(Node) "run()" or \ref init()
    ///must be called before using this function and
    ///node \c v must be reached but not necessarily processed.
    Value currentDist(Node v) const {
      return processed(v) ? (*_dist)[v] : (*_heap)[v];
    }

    ///@}
  };
  ///Default traits class of astar2dir() function.

  ///Default traits class of astar2dir() function.
  ///\tparam GR The type of the digraph.
  ///\tparam LEN The type of the length map.
  template<class GR, class LEN>
  struct Astar2dirWizardDefaultTraits
  {
    ///The type of the digraph the algorithm runs on.
    typedef GR Digraph;
    ///The type of the map that stores the arc lengths.

    ///The type of the map that stores the arc lengths.
    ///It must conform to the \ref concepts::ReadMap "ReadMap" concept.
    typedef LEN LengthMap;
    ///The type of the arc lengths.
    typedef typename LEN::Value Value;

    /// Operation traits for Astar2dir algorithm.

    /// This class defines the operations that are used in the algorithm.
    /// \see Astar2dirDefaultOperationTraits
    typedef Astar2dirDefaultOperationTraits<Value> OperationTraits;

    /// The cross reference type used by the heap.

    /// The cross reference type used by the heap.
    /// Usually it is \c Digraph::NodeMap<int>.
    typedef typename Digraph::template NodeMap<int> HeapCrossRef;
    ///Instantiates a \ref HeapCrossRef.

    ///This function instantiates a \ref HeapCrossRef.
    /// \param g is the digraph, to which we would like to define the
    /// HeapCrossRef.
    static HeapCrossRef *createHeapCrossRef(const Digraph &g)
    {
      return new HeapCrossRef(g);
    }

    ///The heap type used by the Astar2dir algorithm.

    ///The heap type used by the Astar2dir algorithm.
    ///
    ///\sa BinHeap
    ///\sa Astar2dir
    typedef BinHeap<Value, typename Digraph::template NodeMap<int>,
                    std::less<Value> > Heap;

    ///Instantiates a \ref Heap.

    ///This function instantiates a \ref Heap.
    /// \param r is the HeapCrossRef which is used.
    static Heap *createHeap(HeapCrossRef& r)
    {
      return new Heap(r);
    }

    ///\brief The type of the map that stores the predecessor
    ///arcs of the shortest paths.
    ///
    ///The type of the map that stores the predecessor
    ///arcs of the shortest paths.
    ///It must conform to the \ref concepts::WriteMap "WriteMap" concept.
    typedef typename Digraph::template NodeMap<typename Digraph::Arc> PredMap;
    ///Instantiates a PredMap.

    ///This function instantiates a PredMap.
    ///\param g is the digraph, to which we would like to define the
    ///PredMap.
    static PredMap *createPredMap(const Digraph &g)
    {
      return new PredMap(g);
    }

    ///The type of the map that indicates which nodes are processed.

    ///The type of the map that indicates which nodes are processed.
    ///It must conform to the \ref concepts::WriteMap "WriteMap" concept.
    ///By default, it is a NullMap.
    typedef NullMap<typename Digraph::Node,bool> ProcessedMap;
    ///Instantiates a ProcessedMap.

    ///This function instantiates a ProcessedMap.
    ///\param g is the digraph, to which
    ///we would like to define the ProcessedMap.
#ifdef DOXYGEN
    static ProcessedMap *createProcessedMap(const Digraph &g)
#else
    static ProcessedMap *createProcessedMap(const Digraph &)
#endif
    {
      return new ProcessedMap();
    }

    ///The type of the map that stores the distances of the nodes.

    ///The type of the map that stores the distances of the nodes.
    ///It must conform to the \ref concepts::WriteMap "WriteMap" concept.
    typedef typename Digraph::template NodeMap<typename LEN::Value> DistMap;
    ///Instantiates a DistMap.

    ///This function instantiates a DistMap.
    ///\param g is the digraph, to which we would like to define
    ///the DistMap
    static DistMap *createDistMap(const Digraph &g)
    {
      return new DistMap(g);
    }

    ///The type of the shortest paths.

    ///The type of the shortest paths.
    ///It must conform to the \ref concepts::Path "Path" concept.
    typedef lemon::Path<Digraph> Path;
  };

  /// Default traits class used by Astar2dirWizard

  /// Default traits class used by Astar2dirWizard.
  /// \tparam GR The type of the digraph.
  /// \tparam LEN The type of the length map.
  template<typename GR, typename LEN>
  class Astar2dirWizardBase : public Astar2dirWizardDefaultTraits<GR,LEN>
  {
    typedef Astar2dirWizardDefaultTraits<GR,LEN> Base;
  protected:
    //The type of the nodes in the digraph.
    typedef typename Base::Digraph::Node Node;

    //Pointer to the digraph the algorithm runs on.
    void *_g;
    //Pointer to the length map.
    void *_length;
    //Pointer to the map of processed nodes.
    void *_processed;
    //Pointer to the map of predecessors arcs.
    void *_pred;
    //Pointer to the map of distances.
    void *_dist;
    //Pointer to the shortest path to the target node.
    void *_path;
    //Pointer to the distance of the target node.
    void *_di;

  public:
    /// Constructor.

    /// This constructor does not require parameters, therefore it initiates
    /// all of the attributes to \c 0.
    Astar2dirWizardBase() : _g(0), _length(0), _processed(0), _pred(0),
                           _dist(0), _path(0), _di(0) {}

    /// Constructor.

    /// This constructor requires two parameters,
    /// others are initiated to \c 0.
    /// \param g The digraph the algorithm runs on.
    /// \param l The length map.
    Astar2dirWizardBase(const GR &g,const LEN &l) :
      _g(reinterpret_cast<void*>(const_cast<GR*>(&g))),
      _length(reinterpret_cast<void*>(const_cast<LEN*>(&l))),
      _processed(0), _pred(0), _dist(0), _path(0), _di(0) {}

  };

  /// Auxiliary class for the function-type interface of Astar2dir algorithm.

  /// This auxiliary class is created to implement the
  /// \ref astar2dir() "function-type interface" of \ref Astar2dir algorithm.
  /// It does not have own \ref run(Node) "run()" method, it uses the
  /// functions and features of the plain \ref Astar2dir.
  ///
  /// This class should only be used through the \ref astar2dir() function,
  /// which makes it easier to use the algorithm.
  ///
  /// \tparam TR The traits class that defines various types used by the
  /// algorithm.
  template<class TR>
  class Astar2dirWizard : public TR
  {
    typedef TR Base;

    typedef typename TR::Digraph Digraph;

    typedef typename Digraph::Node Node;
    typedef typename Digraph::NodeIt NodeIt;
    typedef typename Digraph::Arc Arc;
    typedef typename Digraph::OutArcIt OutArcIt;

    typedef typename TR::LengthMap LengthMap;
    typedef typename LengthMap::Value Value;
    typedef typename TR::PredMap PredMap;
    typedef typename TR::DistMap DistMap;
    typedef typename TR::ProcessedMap ProcessedMap;
    typedef typename TR::Path Path;
    typedef typename TR::Heap Heap;

  public:

    /// Constructor.
    Astar2dirWizard() : TR() {}

    /// Constructor that requires parameters.

    /// Constructor that requires parameters.
    /// These parameters will be the default values for the traits class.
    /// \param g The digraph the algorithm runs on.
    /// \param l The length map.
    Astar2dirWizard(const Digraph &g, const LengthMap &l) :
      TR(g,l) {}

    ///Copy constructor
    Astar2dirWizard(const TR &b) : TR(b) {}

    ~Astar2dirWizard() {}

    ///Runs Astar2dir algorithm from the given source node.

    ///This method runs %Astar2dir algorithm from the given source node
    ///in order to compute the shortest path to each node.
    void run(Node s)
    {
      Astar2dir<Digraph,LengthMap,TR>
        dijk(*reinterpret_cast<const Digraph*>(Base::_g),
             *reinterpret_cast<const LengthMap*>(Base::_length));
      if (Base::_pred)
        dijk.predMap(*reinterpret_cast<PredMap*>(Base::_pred));
      if (Base::_dist)
        dijk.distMap(*reinterpret_cast<DistMap*>(Base::_dist));
      if (Base::_processed)
        dijk.processedMap(*reinterpret_cast<ProcessedMap*>(Base::_processed));
      dijk.run(s);
    }

    ///Finds the shortest path between \c s and \c t.

    ///This method runs the %Astar2dir algorithm from node \c s
    ///in order to compute the shortest path to node \c t
    ///(it stops searching when \c t is processed).
    ///
    ///\return \c true if \c t is reachable form \c s.
    bool runastar(Node& s, Node& t,ListDigraph::NodeMap<vector<Node>> & _neighbor,
	  ListDigraph::NodeMap<double> & _dists,ListDigraph::NodeMap<double> & _distt,
	   ListDigraph::NodeMap<double> & potential , double shortest_distance, double _INFINITY, ListDigraph::NodeMap<vector<Node>> & _parent, double adaptFactor
 		)
    {
      Astar2dir<Digraph,LengthMap,TR>
        dijk(*reinterpret_cast<const Digraph*>(Base::_g),
             *reinterpret_cast<const LengthMap*>(Base::_length));
      if (Base::_pred)
        dijk.predMap(*reinterpret_cast<PredMap*>(Base::_pred));
      if (Base::_dist)
        dijk.distMap(*reinterpret_cast<DistMap*>(Base::_dist));
      if (Base::_processed)
        dijk.processedMap(*reinterpret_cast<ProcessedMap*>(Base::_processed));
      dijk.runastar(s,t,_neighbor,_dists, _distt,potential, shortest_distance, _INFINITY,_parent, adaptFactor);
      if (Base::_path)
        *reinterpret_cast<Path*>(Base::_path) = dijk.path(t);
      if (Base::_di)
        *reinterpret_cast<Value*>(Base::_di) = dijk.dist(t);
      return dijk.reached(t);
    }
     
bool runastarAdaptive(Node& s, Node& t, ListDigraph::NodeMap<double> & _dists, double heuristicPercent, double _INFINITY)
    {
      Astar2dir<Digraph,LengthMap,TR>
        dijk(*reinterpret_cast<const Digraph*>(Base::_g),
             *reinterpret_cast<const LengthMap*>(Base::_length));
      if (Base::_pred)
        dijk.predMap(*reinterpret_cast<PredMap*>(Base::_pred));
      if (Base::_dist)
        dijk.distMap(*reinterpret_cast<DistMap*>(Base::_dist));
      if (Base::_processed)
        dijk.processedMap(*reinterpret_cast<ProcessedMap*>(Base::_processed));
      dijk.runastarAdaptive(s,t,_dists, heuristicPercent, _INFINITY);
      if (Base::_path)
        *reinterpret_cast<Path*>(Base::_path) = dijk.path(t);
      if (Base::_di)
        *reinterpret_cast<Value*>(Base::_di) = dijk.dist(t);
      return dijk.reached(t);
    }

    template<class T>
    struct SetPredMapBase : public Base {
      typedef T PredMap;
      static PredMap *createPredMap(const Digraph &) { return 0; };
      SetPredMapBase(const TR &b) : TR(b) {}
    };

    ///\brief \ref named-templ-param "Named parameter" for setting
    ///the predecessor map.
    ///
    ///\ref named-templ-param "Named parameter" function for setting
    ///the map that stores the predecessor arcs of the nodes.
    template<class T>
    Astar2dirWizard<SetPredMapBase<T> > predMap(const T &t)
    {
      Base::_pred=reinterpret_cast<void*>(const_cast<T*>(&t));
      return Astar2dirWizard<SetPredMapBase<T> >(*this);
    }

    template<class T>
    struct SetDistMapBase : public Base {
      typedef T DistMap;
      static DistMap *createDistMap(const Digraph &) { return 0; };
      SetDistMapBase(const TR &b) : TR(b) {}
    };

    ///\brief \ref named-templ-param "Named parameter" for setting
    ///the distance map.
    ///
    ///\ref named-templ-param "Named parameter" function for setting
    ///the map that stores the distances of the nodes calculated
    ///by the algorithm.
    template<class T>
    Astar2dirWizard<SetDistMapBase<T> > distMap(const T &t)
    {
      Base::_dist=reinterpret_cast<void*>(const_cast<T*>(&t));
      return Astar2dirWizard<SetDistMapBase<T> >(*this);
    }

    template<class T>
    struct SetProcessedMapBase : public Base {
      typedef T ProcessedMap;
      static ProcessedMap *createProcessedMap(const Digraph &) { return 0; };
      SetProcessedMapBase(const TR &b) : TR(b) {}
    };

    ///\brief \ref named-func-param "Named parameter" for setting
    ///the processed map.
    ///
    ///\ref named-templ-param "Named parameter" function for setting
    ///the map that indicates which nodes are processed.
    template<class T>
    Astar2dirWizard<SetProcessedMapBase<T> > processedMap(const T &t)
    {
      Base::_processed=reinterpret_cast<void*>(const_cast<T*>(&t));
      return Astar2dirWizard<SetProcessedMapBase<T> >(*this);
    }

    template<class T>
    struct SetPathBase : public Base {
      typedef T Path;
      SetPathBase(const TR &b) : TR(b) {}
    };

    ///\brief \ref named-func-param "Named parameter"
    ///for getting the shortest path to the target node.
    ///
    ///\ref named-func-param "Named parameter"
    ///for getting the shortest path to the target node.
    template<class T>
    Astar2dirWizard<SetPathBase<T> > path(const T &t)
    {
      Base::_path=reinterpret_cast<void*>(const_cast<T*>(&t));
      return Astar2dirWizard<SetPathBase<T> >(*this);
    }

    ///\brief \ref named-func-param "Named parameter"
    ///for getting the distance of the target node.
    ///
    ///\ref named-func-param "Named parameter"
    ///for getting the distance of the target node.
    Astar2dirWizard dist(const Value &d)
    {
      Base::_di=reinterpret_cast<void*>(const_cast<Value*>(&d));
      return *this;
    }

  };

  
  ///Function-type interface for Astar2dir algorithm.

  /// \ingroup shortest_path
  ///Function-type interface for Astar2dir algorithm.
  ///
  ///This function also has several \ref named-func-param "named parameters",
  ///they are declared as the members of class \ref Astar2dirWizard.
  ///The following examples show how to use these parameters.
  ///\code
  ///  // Compute shortest path from node s to each node
  ///  astar2dir(g,length).predMap(preds).distMap(dists).run(s);
  ///
  ///  // Compute shortest path from s to t
  ///  bool reached = astar2dir(g,length).path(p).dist(d).run(s,t);
  ///\endcode
  ///\warning Don't forget to put the \ref Astar2dirWizard::run(Node) "run()"
  ///to the end of the parameter list.
  ///\sa Astar2dirWizard
  ///\sa Astar2dir
  template<typename GR, typename LEN>
  Astar2dirWizard<Astar2dirWizardBase<GR,LEN> >
  astar2dir(const GR &digraph, const LEN &length)
  {
    return Astar2dirWizard<Astar2dirWizardBase<GR,LEN> >(digraph,length);
  }
  

} //END OF NAMESPACE LEMON

#endif
