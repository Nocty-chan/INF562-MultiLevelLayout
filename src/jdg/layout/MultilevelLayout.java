package jdg.layout;

import java.util.ArrayList;
import java.util.HashSet;
import java.util.LinkedList;
import java.util.List;
import java.util.Queue;
import java.util.Random;
import java.util.Set;
import java.util.Stack;

import Jcg.geometry.Point_3;
import Jcg.geometry.Vector_3;
import jdg.graph.AdjacencyListGraph;
import jdg.graph.Node;

/**
 * A class implementing the multi-level force directed algorithm by Walshaw (2003)
 * 
 * @author Luca Castelli Aleardi, Ecole Polytechnique
 * @version fev 2017
 */
public class MultilevelLayout extends FR91Layout {
	// parameters of the algorithm
	public static int threshold = 10; //minimal size for coarsening process
	public LinkedList<AdjacencyListGraph> graphs; // sequence of coarser graphs
	private static Random randomInt = new Random();
	/**
	 * Initialize the parameters of the force-directed layout
	 * 
	 *  @param g  input graph to draw
	 *  @param w  width of the drawing area
	 *  @param h  height of the drawing area
	 */
	public MultilevelLayout(AdjacencyListGraph g, double w, double h) {
	  super(g, w, h);
		System.out.print("Initializing force-directed method (Walshaw 2003)...");
		if(g==null) {
			System.out.println("Input graph not defined");
			System.exit(0);
		}
		this.g=g;
		int N=g.sizeVertices();
		
		// set the parameters of the algorithm
		
		System.out.println("done ("+N+" nodes)");
		System.out.println("k="+k+" - temperature="+temperature);
		System.out.println(this.toString());
	}
	
	public AdjacencyListGraph getGraph(int selectedGraph) {
	  if (this.graphs!= null && selectedGraph < this.graphs.size()) {
	    return this.graphs.get(selectedGraph);
	  }
	  return this.g;
	}
	/**
	 * Enable cooling process
	 */	
	public void enableCooling() {
		this.useCooling=true;
	}

	/**
	 * Enable cooling process
	 */	
	public void disableCooling() {
		this.useCooling=false;
	}

	/**
	 * Compute the (intensity of the) attractive force between two nodes at a given distance
	 * 
	 * @param distance  distance between two nodes
	 */	
	public double attractiveForce(double distance) {
		return (distance*distance)/k;
	}
	
	/**
	 * Compute the (intensity of the) repulsive force between two nodes at a given distance
	 * 
	 * @param distance  distance between two nodes
	 */	
	public double repulsiveForce(double distance, double weight) {
		return C * weight *  (k*k)/distance;
	}

	/**
	 * Perform the multi-level Force-Directed algorithm.
	 * Positions of vertices are updated according to their mutual attractive and repulsive forces.
	 */	
	public void computeLayout(int selectedGraph) {
		if(iterationCount>=maxIterations)
			return;
		if (selectedGraph < this.graphs.size()) {
		  System.out.println(selectedGraph);
	    this.computeLayoutOneGraph(this.graphs.get(selectedGraph));  
		} else {
		  this.computeLayoutOneGraph(this.graphs.get(this.graphs.size() - 1));
		}
		
	}

	public void simplify() {
	  coarsenGraph();
	}
	
	 public void coarsenGraph() {
	   this.graphs = new LinkedList<AdjacencyListGraph>();
	   this.graphs.add(this.g);
	   // Set weights */
	   for (Node v : this.g.vertices) {
	     v.weight = 1.0d;
	   }
	   AdjacencyListGraph coarserGraph = simplify(this.g);
	   while (coarserGraph != null) {
	     this.graphs.add(coarserGraph);
	     coarserGraph = simplify(coarserGraph);
	   }
	   System.out.println("Coarsening process created "+ this.graphs.size() + "graphs.");
	 }
	/**
	 * Given g, compute a coarser graph g', obtained performing a maximal sequence of edge collapses
	 */	
	public static AdjacencyListGraph simplify(AdjacencyListGraph g) {
	  System.out.println("Simplifying graph");
	  if (g.sizeVertices() < threshold) return null;
		AdjacencyListGraph coarserGraph = g.getCopy();
    /* Initialize all nodes
     * By default, tag is 0 and descendant is itself in new graph.
     */
		Stack<Node> nodes = new Stack<Node>();
		for (Node v: coarserGraph.vertices) {
		  nodes.add(v);
		  v.tag = 0;
		  g.getNode(v.getLabel()).descendant = v;
		}
		
		int count = 0;
		/* Iterate over all nodes of coarserGraph */
		while(!nodes.isEmpty()) {
		  int index = randomInt.nextInt(nodes.size());
		  Node v = nodes.get(index);
		  nodes.remove(index);
		  if (v.tag == 1) continue;
		  count++;
		  /* Get minimal weight neighbour that is unmarked */
		  double minWeight = Double.MAX_VALUE;
		  Node minWeightNeighbour = null;
		  for (Node u : v.neighborsList()) {
		    if (u.tag == 0 && minWeight > u.weight) {
		      minWeight = u.weight;
		      minWeightNeighbour = u;
		    }
		  }
		  if (minWeightNeighbour != null) {
		    //Collapse (minWeightNeighbour,v).

	      System.out.println("Collapsing edge ("+minWeightNeighbour.getLabel() + v.getLabel());
		    //Update neighbours list.
		    ArrayList<Node> neighbours = (ArrayList<Node>) minWeightNeighbour.neighbors.clone();
		    for (Node neighbour : neighbours) {
		      if (!neighbour.equals(v)) {
		        coarserGraph.addEdge(v, neighbour);
		        coarserGraph.addEdge(neighbour, v);
	        }
		    }
		    
        //Remove node 
        coarserGraph.removeNode(minWeightNeighbour);
  
        // Update tag for all neighbours of v

        for (Node neighbour : v.neighbors) {
		      neighbour.tag = 1;
		    }
		    
		    // Update descendant and tag for collapsed vertices.
        g.getNode(minWeightNeighbour.getLabel()).descendant = v;
        minWeightNeighbour.tag = 1;
		    v.tag = 1;
		    
		    //update weights
		    v.weight +=minWeightNeighbour.weight;
		  }
		}
	return coarserGraph;	
	}

	public String toString() {
		String result="Multi-level force-directed algorihm (Walshaw 2003)\n";

		return result;
	}
	
	
	public void computeLayoutOneGraph(AdjacencyListGraph graph) {
    if(iterationCount>=maxIterations)
      return;
    
    System.out.print("Performing iteration: "+this.iterationCount);
    if(useCooling==true)
      System.out.println("(temperature: "+this.temperature+")");
    else
      System.out.println("(no simulated annealing)");
        
    long startTime=System.nanoTime(), endTime; // for evaluating time performances
    
    // first step: for each vertex compute the displacements due to attractive and repulsive forces
    Vector_3[] attractiveDisp=this.computeAllAttractiveForces(graph); // displacement due to attractive forces between neighbors
    Vector_3[] repulsiveDisp=this.computeAllRepulsiveForces(graph); // displacement due to repulsive forces between all nodes
    
    // second step: compute the total displacements and move all nodes to their new locations
    for(int i=0;i<graph.sizeVertices();i++) {
      Node u=graph.vertices.get(i);
      if(u.degree()>0) { 
        // only non isolated nodes are involved in the computation
        Point_3 p=u.getPoint(); // geometric position of node 'u'

        Vector_3 disp=attractiveDisp[i].sum(repulsiveDisp[i]); // total displacement for vertex 'u' 
        double norm=Math.sqrt(disp.squaredLength().doubleValue()); // length of the displacement of vertex 'u'
        
        // use temperature to limit the displacement
        double min=Math.min(norm, this.temperature);
        disp=disp.multiplyByScalar(min/norm);
        
        p.translateOf(disp);
      }
    }
    
    // evaluate time performances
      endTime=System.nanoTime();
        double duration=(double)(endTime-startTime)/1000000000.;
        System.out.println("iteration "+this.iterationCount+" done ("+duration+" seconds)");

    if(useCooling==true) // check whether "simulated annealing" is used
      this.cooling(); // update temperature
    
    this.iterationCount++; // increase counter (to count the number of performed iterations)
  }
  
  /**
   * Compute the displacement of vertex 'u', due to repulsive forces (of all nodes)
   * 
   * @param u  the vertex to which repulsive forces are applied
   * @return 'displacement' a 3d vector storing the displacement of vertex 'u'
   */ 
  public Vector_3 computeRepulsiveForce(Node u, AdjacencyListGraph graph) {
    Vector_3 displacement=new Vector_3(0., 0., 0.);
    Vector_3 delta;
    Point_3 p=u.getPoint();
    
    for(Node v: graph.vertices) {
      if(v!=null && u!=v && v.degree()>=0) {
        delta=new Vector_3(v.getPoint(),p);
        double norm=Math.sqrt(delta.squaredLength().doubleValue());
        displacement=displacement.sum(delta.multiplyByScalar(repulsiveForce(norm, v.weight)));
      }
    }

    return displacement;
  }
  
  /**
   * Compute, for each vertex, the displacement due to repulsive forces (between all nodes)
   * 
   * @return a vector v[]: v[i] stores the geometric displacement of the i-th node
   */ 
  public Vector_3[] computeAllRepulsiveForces(AdjacencyListGraph graph) {
    Vector_3[] repulsiveDisp=new Vector_3[graph.sizeVertices()];
    for(int i=0;i<graph.vertices.size();i++)
      repulsiveDisp[i]=new Vector_3(0., 0., 0.);
    
    long startTime=System.nanoTime(), endTime; // for evaluating time performances
    
    int countVertices=0;
    for(Node u: graph.vertices) {
      //repulsiveForce(u, repulsiveDisp); // faster version (forces are computed only once)
      repulsiveDisp[countVertices]=computeRepulsiveForce(u, graph); // slow version (forces are computed twice)
      countVertices++;
    }
      endTime=System.nanoTime();
        double duration=(double)(endTime-startTime)/1000000000.;
        System.out.print("\ttimings: "+duration+" seconds (repulsive force)");
        System.out.println(" ["+this.countRepulsive+" forces computed]");
        this.countRepulsive=0;
        
        return repulsiveDisp;
  }
  
  /**
   * Compute the displacement of vertex 'u', due to the attractive forces of its neighbors
   * 
   * @param u  the vertex to which attractive forces are applied
   * @return 'disp' a 3d vector storing the displacement of vertex 'u'
   */ 
  public Vector_3 computeAttractiveForce(Node u) {
    Vector_3 displacement=new Vector_3(0., 0., 0.);
    Vector_3 delta;
    Point_3 p=u.getPoint();
    for(Node v: u.neighbors) {
      if(v!=null) {
        delta=new Vector_3(p, v.getPoint());
        double norm=Math.sqrt(delta.squaredLength().doubleValue());
        displacement=displacement.sum(delta.multiplyByScalar(attractiveForce(norm)));
      }
    }
    return displacement;
  }
  
  /**
   * Compute, for each vertex, the displacement due to attractive forces (between neighboring nodes)
   * 
   * @return a vector v[]: v[i] stores the geometric displacement of the i-th node
   */ 
  public Vector_3[] computeAllAttractiveForces(AdjacencyListGraph graph) {
    Vector_3[] attractiveDisp=new Vector_3[graph.sizeVertices()];
    long startTime=System.nanoTime(), endTime; // for evaluating time performances
    
    int i=0;
    for(Node u: graph.vertices) {
      attractiveDisp[i]=computeAttractiveForce(u);
      i++;
    }
      endTime=System.nanoTime();
        double duration=(double)(endTime-startTime)/1000000000.;
        System.out.println("\ttimings: "+duration+" seconds (attractive force)");
        
        return attractiveDisp;
  }
  
}
