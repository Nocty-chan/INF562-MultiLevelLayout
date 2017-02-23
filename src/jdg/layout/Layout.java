package jdg.layout;

import java.util.Random;

import jdg.graph.AdjacencyListGraph;
import jdg.graph.Node;

import Jcg.geometry.Point_2;
import Jcg.geometry.Point_3;
import Jcg.geometry.Transformation_2;
import Jcg.geometry.Vector_2;
import Jcg.viewer.old.Fenetre;

/**
 * Network visualization: abstract implementation of the force-directed paradigm
 * 
 * @author Luca Castelli Aleardi, Ecole Polytechnique
 * @version august 2015
 */
public abstract class Layout {
	public AdjacencyListGraph g;
	public double w, h; // dimensions of the drawing area
	
	public static double coolingConstant=0.50;

	public static int maxIterations=1000; // maximum number of iterations to perform
    
	public static int seed=10;
	/** Random generator */	
	static Random generator = new Random(seed); // initialize random generator
	
	/**
	 * Initialize vertex locations at random (in a square of given size WxH)
	 */	
	public static void setRandomPoints(AdjacencyListGraph g, double w, double h) {
		Point_3 p;
		double w1=w/2., h1=h/2.;
		for(Node u: g.vertices){
			double n1=w1-2*w1*generator.nextDouble();
			double n2=h1-2*h1*generator.nextDouble();
		    p = new Point_3 (n1, n2, 0.0);
		    u.setPoint(p);
		}		
	}

	/**
	 * Perform one iteration of the Force-Directed algorithm.
	 * Positions of vertices are updated according to their
	 * mutual attractive and repulsive forces.
	 */	
	public abstract void computeLayout();
	
	/**
	 * Enable cooling process
	 */	
	public abstract void enableCooling();
	
	/**
	 * Enable cooling process
	 */	
	public abstract void disableCooling();
		
	
	public abstract void simplify();
	/**
	 * draw the graph in a 2D frame
	 */
	public void draw2D() {
		Fenetre f=new Fenetre();
		for(Node u: this.g.vertices) {
			Point_2 p=new Point_2(u.getPoint().getX(), u.getPoint().getY());
			for(Node v: u.neighbors) {
				if(v!=null && v.index>u.index) { // draw only directed edges (u, v) such that u<v
					Point_2 q=new Point_2(v.getPoint().getX(), v.getPoint().getY());
					f.addSegment(p, q);
				}
			}
		}
	}
	
}
