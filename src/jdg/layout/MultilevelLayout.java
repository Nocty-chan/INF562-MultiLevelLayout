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
	public static int threshold; //minimal size for coarsening process
	public LinkedList<AdjacencyListGraph> graphs; // sequence of coarser graphs
	private static Random randomInt = new Random();
	private static int totalNumberOfIterations;
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
		this.threshold = this.g.sizeVertices()/5;
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
	public double attractiveForce(double distance, int index) {
	  double updatedConstant = k * Math.pow(Math.sqrt(7.0/4.0), index);
    return (distance*distance)/updatedConstant;
	}
	
	/**
	 * Compute the (intensity of the) repulsive force between two nodes at a given distance
	 * 
	 * @param distance  distance between two nodes
	 */	
	public double repulsiveForce(double distance, double weight, int index) {
	  double updatedConstant = k * Math.pow(Math.sqrt(7.0/4.0), index);
		return C * weight * (updatedConstant * updatedConstant)/distance;
	}

	/**
	 * Perform the multi-level Force-Directed algorithm.
	 * Positions of vertices are updated according to their mutual attractive and repulsive forces.
	 */	
	public void computeLayout(int selectedGraph) {
	  if (this.graphs == null) 
	    coarsenGraph();
		if(iterationCount>=maxIterations)
			return;
		for (int i = this.graphs.size() - 1; i > 0; i--) {
		  this.temperature = this.w / 10;
		  this.computeLayoutOneGraph(i, this.graphs.get(i));
		  //Update position
		  for (Node u: this.graphs.get(i - 1).vertices) {
		   Point_3 p = u.p;
		    Point_3 descendant = (u.descendant == null) ? (this.graphs.get(i).getNode(u.getLabel()).p): u.descendant.p;
		    p.x = descendant.x.doubleValue();
        p.y = descendant.y.doubleValue();
        p.z = descendant.z.doubleValue();
		    if (u.descendant != null) {
		      u.p.x += 10 * Math.min(1.0, randomInt.nextInt(20));
		      u.p.y += 10 * Math.min(1.0, randomInt.nextInt(20));
		      u.p.z += 10 * Math.min(1.0, randomInt.nextInt(20));
		    }		      
		  }		  
		}
		this.computeLayoutOneGraph(this.graphs.get(0));
		this.iterationCount++;
		System.out.println("Found layout in "+ this.totalNumberOfIterations + " iterations in total.");
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
	
	
	public void computeLayoutOneGraph(int index, AdjacencyListGraph graph) {
    if(useCooling==true)
      System.out.println("(temperature: "+this.temperature+")");
    else
      System.out.println("(no simulated annealing)");
    int converged = 0;
    while (converged != 1) {
      converged = 1;
      for (Node u : graph.vertices) {
        Vector_3 repulsiveDisplacement = this.computeRepulsiveForce(u, graph, index);
        Vector_3 attractiveDisplacement = this.computeAttractiveForce(u, index);
        Vector_3 displacement = repulsiveDisplacement.sum(attractiveDisplacement);
        double norm = Math.sqrt(displacement.squaredLength().doubleValue());
        double min = Math.min(norm, this.temperature);
        displacement = displacement.multiplyByScalar(min/norm);
        u.getPoint().translateOf(displacement);
        
        if (displacement.squaredLength().doubleValue() > k) {
          converged = 0;
        }
      }
      if(this.useCooling) {
        this.cooling();
      }
     this.totalNumberOfIterations++;
    }
    System.out.println("Computed layout of graph "+ index);
  }
  
  /**
   * Compute the displacement of vertex 'u', due to repulsive forces (of all nodes)
   * 
   * @param u  the vertex to which repulsive forces are applied
   * @return 'displacement' a 3d vector storing the displacement of vertex 'u'
   */ 
  public Vector_3 computeRepulsiveForce(Node u, AdjacencyListGraph graph, int index) {
    Vector_3 displacement=new Vector_3(0., 0., 0.);
    Vector_3 delta;
    Point_3 p=u.getPoint();
    
    for(Node v: graph.vertices) {
      if(v!=null && u!=v && v.degree()>=0) {
        delta=new Vector_3(v.getPoint(),p);
        double norm=Math.sqrt(delta.squaredLength().doubleValue());
        displacement=displacement.sum(delta.multiplyByScalar(repulsiveForce(norm, v.weight, index)));
      }
    }

    return displacement;
  }
  
  /**
   * Compute the displacement of vertex 'u', due to the attractive forces of its neighbors
   * 
   * @param u  the vertex to which attractive forces are applied
   * @return 'disp' a 3d vector storing the displacement of vertex 'u'
   */ 
  public Vector_3 computeAttractiveForce(Node u, int index) {
    Vector_3 displacement=new Vector_3(0., 0., 0.);
    Vector_3 delta;
    Point_3 p=u.getPoint();
    for(Node v: u.neighbors) {
      if(v!=null) {
        delta=new Vector_3(p, v.getPoint());
        double norm=Math.sqrt(delta.squaredLength().doubleValue());
        displacement=displacement.sum(delta.multiplyByScalar(attractiveForce(norm, index)));
      }
    }
    return displacement;
  }  
}
