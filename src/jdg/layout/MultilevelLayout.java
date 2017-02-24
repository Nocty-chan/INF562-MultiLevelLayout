package jdg.layout;

import java.util.ArrayList;
import java.util.HashSet;
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
public class MultilevelLayout extends Layout {
	// parameters of the algorithm
	public int iterationCount=0; // count the number of performed iterations
	public double k; // natural spring length
	public double area; // area of the drawing (width times height)
	public double C; // step
	public double gamma; //
	public boolean useCooling; // say whether performing simulated annealing
	public static int threshold = 2; //minimal size for coarsening process
	public AdjacencyListGraph[] graphs; // sequence of coarser graphs
	private static Random randomInt = new Random();

	/**
	 * Initialize the parameters of the force-directed layout
	 * 
	 *  @param g  input graph to draw
	 *  @param w  width of the drawing area
	 *  @param h  height of the drawing area
	 */
	public MultilevelLayout(AdjacencyListGraph g, double w, double h) {
		System.out.print("Initializing force-directed method (Walshaw 2003)...");
		if(g==null) {
			System.out.println("Input graph not defined");
			System.exit(0);
		}
		this.g=g;
		int N=g.sizeVertices();
		
		// set the parameters of the algorithm
		
		System.out.println("done ("+N+" nodes)");
		//System.out.println("k="+k+" - temperature="+temperature);
		System.out.println(this.toString());
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
		return -C * weight *  (k*k)/distance;
	}

	/**
	 * Perform the multi-level Force-Directed algorithm.
	 * Positions of vertices are updated according to their mutual attractive and repulsive forces.
	 */	
	public void computeLayout() {
		if(iterationCount>=maxIterations)
			return;

		throw new Error("To be completed");		
	}
	
	public void simplify() {
	  
	}

	/**
	 * Given g, compute a coarser graph g', obtained performing a maximal sequence of edge collapses
	 */	
	public static AdjacencyListGraph simplify(AdjacencyListGraph g) {
	  System.out.println("Simplifying graph");
	  if (g.sizeVertices() < threshold) return g;
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
	
}
