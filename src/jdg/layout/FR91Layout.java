package jdg.layout;

import jdg.graph.AdjacencyListGraph;
import jdg.graph.Node;

import Jcg.geometry.Point_3;
import Jcg.geometry.Vector_3;

/**
 * A class implementing the force directed algorithm by Fruchterman and Reingold (1991)
 * 
 * @author Luca Castelli Aleardi, Ecole Polytechnique
 * @version fev 2017
 */
public class FR91Layout extends Layout {
	// parameters of the algorithm by Fruchterman and Reingold
	public double k; // natural spring length
	public double area; // area of the drawing (width times height)
	public double C; // step
	public double temperature; // initial temperature
	public double minTemperature; // minimal temperature (strictly positive)
	public double coolingConstant; // constant term: the temperature decreases linearly at each iteration
	public boolean useCooling; // say whether performing simulated annealing
	
	public int iterationCount=0; // count the number of performed iterations
	private int countRepulsive=0; // count the number of computed repulsive forces (to measure time performances)
	
	/**
	 * Initialize the parameters of the force-directed layout
	 * 
	 *  @param g  input graph to draw
	 *  @param w  width of the drawing area
	 *  @param h  height of the drawing area
	 *  @param C  step length
	 */
	public FR91Layout(AdjacencyListGraph g, double w, double h) {
		System.out.print("Initializing force-directed method: Fruchterman-Reingold 91...");
		if(g==null) {
			System.out.println("Input graph not defined");
			System.exit(0);
		}
		this.g=g;
		int N=g.sizeVertices();
		
		// set the parameters of the algorithm FR91
		this.C=1.;
		this.w=w;
		this.h=h;
		this.area=w*h;
		this.k=C*Math.sqrt(area/N);
		this.temperature=w/10.; // the temperature is a fraction of the width of the drawing area
		this.minTemperature=0.01;
		this.coolingConstant=0.50;
		this.useCooling=true; // use cooling system by default
		
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
	public double repulsiveForce(double distance) {
		countRepulsive++;
		return (k*k)/distance;
	}

	public void simplify() {
	  AdjacencyListGraph newGraph = MultilevelLayout.simplify(this.g);
	  this.g.vertices = newGraph.vertices;
	}
	/**
	 * Perform one iteration of the Force-Directed algorithm.
	 * Positions of vertices are updated according to their mutual attractive and repulsive forces.
	 */	
	public void computeLayout() {
		if(iterationCount>=maxIterations)
			return;
		
		System.out.print("Performing iteration: "+this.iterationCount);
		if(useCooling==true)
			System.out.println("(temperature: "+this.temperature+")");
		else
			System.out.println("(no simulated annealing)");
				
		long startTime=System.nanoTime(), endTime; // for evaluating time performances
		
		// first step: for each vertex compute the displacements due to attractive and repulsive forces
		Vector_3[] attractiveDisp=this.computeAllAttractiveForces(); // displacement due to attractive forces between neighbors
		Vector_3[] repulsiveDisp=this.computeAllRepulsiveForces(); // displacement due to repulsive forces between all nodes
		
		// second step: compute the total displacements and move all nodes to their new locations
		for(int i=0;i<this.g.sizeVertices();i++) {
			Node u=this.g.vertices.get(i);
			if(u.degree()>0) { 
			  System.out.println("Non isolated node: " + i + ".");
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
	public Vector_3 computeRepulsiveForce(Node u) {
		Vector_3 displacement=new Vector_3(0., 0., 0.);
		Vector_3 delta;
		Point_3 p=u.getPoint();
		
		for(Node v: this.g.vertices) {
			if(v!=null && u!=v && v.degree()>=0) {
				delta=new Vector_3(v.getPoint(),p);
				double norm=Math.sqrt(delta.squaredLength().doubleValue());
				displacement=displacement.sum(delta.multiplyByScalar(repulsiveForce(norm)));
			}
		}

		return displacement;
	}
	
	/*public void repulsiveForce(Node u, Vector_3[] disp) {
		double px=u.getPoint().x;
		double py=u.getPoint().y;
		double pz=u.getPoint().z;
		
		int indexU=u.index;
		for(int i=0;i<indexU;i++) {
			Node v=this.g.vertices.get(i);
			if(v!=null && v.degree()>0) {
				//delta=new Vector_3(v.getPoint(),p);
				double dx=px-v.getPoint().x;
				double dy=py-v.getPoint().y;
				double dz=pz-v.getPoint().z;
				//double norm=Math.sqrt(delta.squaredLength().doubleValue());
				double norm=Math.sqrt(dx*dx+dy*dy+dz*dz);
				double force=repulsiveForce(norm);
				
				dx=dx*force;
				dy=dy*force;
				dz=dz*force;
				//Vector_3 delta2=delta.multiplyByScalar(force);
				//disp[indexU]=disp[indexU].sum(delta2);
				disp[indexU].x=disp[indexU].x+dx;
				disp[indexU].y=disp[indexU].y+dy;
				disp[indexU].z=disp[indexU].z+dz;

				disp[i].x=disp[i].x+(-1*dx);
				disp[i].y=disp[i].y+(-1*dy);
				disp[i].z=disp[i].z+(-1*dz);
			}
		}
	}*/
	
	/**
	 * Compute, for each vertex, the displacement due to repulsive forces (between all nodes)
	 * 
	 * @return a vector v[]: v[i] stores the geometric displacement of the i-th node
	 */	
	public Vector_3[] computeAllRepulsiveForces() {
		Vector_3[] repulsiveDisp=new Vector_3[this.g.sizeVertices()];
		for(int i=0;i<this.g.vertices.size();i++)
			repulsiveDisp[i]=new Vector_3(0., 0., 0.);
		
		long startTime=System.nanoTime(), endTime; // for evaluating time performances
		
		int countVertices=0;
		for(Node u: this.g.vertices) {
			//repulsiveForce(u, repulsiveDisp); // faster version (forces are computed only once)
			repulsiveDisp[countVertices]=computeRepulsiveForce(u); // slow version (forces are computed twice)
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
	public Vector_3[] computeAllAttractiveForces() {
		Vector_3[] attractiveDisp=new Vector_3[this.g.sizeVertices()];
		long startTime=System.nanoTime(), endTime; // for evaluating time performances
		
		int i=0;
		for(Node u: this.g.vertices) {
			attractiveDisp[i]=computeAttractiveForce(u);
			i++;
		}
    	endTime=System.nanoTime();
        double duration=(double)(endTime-startTime)/1000000000.;
        System.out.println("\ttimings: "+duration+" seconds (attractive force)");
        
        return attractiveDisp;
	}
	
	/**
	 * Cooling system: the temperature decreases linearly at each iteration
	 * 
	 * Remark: the temperature is assumed to remain strictly positive (>=minTemperature)
	 */	
	protected void cooling() {
		this.temperature=Math.max(this.temperature-coolingConstant, minTemperature);
		//this.temperature=this.temperature*0.98; // variant
	}
	
	public String toString() {
		String result="force-directed algorihm: Fruchterman Reingold\n";
		result=result+"\t area= "+w+" x "+h+"\n";
		result=result+"\t k= "+this.k+"\n";
		result=result+"\t C= "+this.C+"\n";
		result=result+"\t initial temperature= "+this.temperature+"\n";
		result=result+"\t minimal temperature= "+this.minTemperature+"\n";
		result=result+"\t cooling constant= "+this.coolingConstant+"\n";
		
		return result;
	}
	
}
