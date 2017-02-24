package jdg;

import java.awt.BorderLayout;

import processing.core.*;

import jdg.graph.AdjacencyListGraph;
import jdg.graph.Node;
import jdg.layout.*;

import Jcg.geometry.Point_2;
import Jcg.geometry.Point_3;
import Jcg.geometry.Vector_2;

/**
 * A class for drawing (dynamic) graphs (using Processing 1.5.1)
 *
 * @author Luca Castelli Aleardi (Ecole Polytechnique, 2015)
 */
public class DrawGraph extends PApplet {
    // coordinates of the bounding box
    protected double xmin=Double.MAX_VALUE, xmax=Double.MIN_VALUE, ymin=Double.MAX_VALUE, ymax=Double.MIN_VALUE;
    private int selectedLayout = 0;
    // parameters for edge rendering
    double boundaryThickness=0.5;
    private int backgroundColor=255;
    private int edgeColor=50;
    private int edgeOpacity=200;
        
    /** node selected with mouse click (to show)  */
    public Node selectedNode=null; 
	  public Point_2 current; // coordinates of the selected point
    
    /** Layout algorithm  */
    public Layout layout;
    
    /** input graph to draw */
    static public AdjacencyListGraph inputGraph;
        	
   	// parameters of the 2d frame/canvas
    public static int sizeX=600; // horizontal size of the canvas (pizels)
    public static int sizeY=600; // vertical size of the canvas (pixels)
    Point_2 a,b; // range of the window (left bottom and right top corners)
	  
	  /**
	   * Initialize the frame
	   */
	  public void setup() {
		  System.out.println("Setting Canvas size: "+sizeX+" x "+sizeY);
		  this.size(sizeX,sizeY); // set the size of the Java Processing frame
		  
		  // set drawing parameters (size and range of the drawing layout)
		  double w2=sizeX/2.0;
		  double h2=sizeY/2.0;
		  this.a=new Point_2(-w2, -h2); // left bottom corner (the drawing region is centered at the origin)
		  this.b=new Point_2(w2, h2); // top right corner of the drawing region
		    	      
	      // set the graph layout method
	    //this.layout=new FR91Layout(inputGraph, sizeX, sizeY); // force-directed method (Fruchterman Reingold)
		  this.layout = new MultilevelLayout(inputGraph, sizeX, sizeY);
		  // at the beginning the nodes are placed at random
		  this.layout.setRandomPoints(inputGraph, w2, h2);
	  }

	  /**
	   * Main method for drawing the applet
	   */
	  public void draw() {
	    this.background(this.backgroundColor); // set the color of background
	    
	    this.display2D(); // draw all edges in gray

	    if(this.selectedNode!=null) {
	    	this.drawVertexLabel(this.selectedNode);
	    }
	    
	    this.drawOptions();
	  }
	  
	  /**
	   * Deal with keyboard events
	   */
	  public void keyPressed(){
		  switch(key) {
		  	case('z'):this.updateBoundingBox(); break;
		  	case('c'):this.layout.computeLayout(); break;
		  	case('o'):this.zoom(1.2); break;
		  	case('i'):this.zoom(0.8); break;
		  	case('l'):this.layout.draw2D(); break;
		  	case('e'):this.layout.enableCooling();; break;
		  	case('d'):this.layout.disableCooling();; break;
		  	case('p'):this.layout.simplify();;break;
		  	case('r'):
		  	  this.selectedLayout++;
          this.inputGraph = this.layout.getGraph(selectedLayout);;break;
		  	case('t'):
		  	  this.selectedLayout= Math.max(0, this.selectedLayout -1);
		  	   this.inputGraph = this.layout.getGraph(selectedLayout); break;
		  }
	  }
	  
	  public void zoom(double factor) {
		  Point_2 barycenter=Point_2.midPoint(a, b);
		  Vector_2 vA=(Vector_2)barycenter.minus(a);
		  Vector_2 vB=(Vector_2)barycenter.minus(b);
		  vA=vA.multiplyByScalar(factor);
		  vB=vB.multiplyByScalar(factor);
		  a=barycenter.sum(vA);
		  b=barycenter.sum(vB);
	  }
	  	  
	  public void mouseClicked() {
		  if(mouseButton==LEFT) { // select a vertex (given its 2D position)
			  this.selectedNode=this.selectNode(mouseX, mouseY);
			  if(selectedNode!=null)
				  System.out.println("vertex "+selectedNode.index);
		  }
	  }

	  public void mousePressed() {
		  this.current=new Point_2(mouseX, mouseY);
	  }
	  
	  public void mouseReleased() {
	  }
	  
	  public void mouseDragged() {
		  if(mouseButton==RIGHT) { // translate the window
			  double norm=Math.sqrt(this.a.squareDistance(this.b).doubleValue());
			  double scaleFactor=norm/(this.sizeX);
			  
			  double deltaX=(mouseX-current.getX().doubleValue())*(scaleFactor);
			  double deltaY=(current.getY().doubleValue()-mouseY)*(scaleFactor);
		  
			  this.a.translateOf(new Vector_2(-deltaX, -deltaY)); // update the left bottom and right top vertices
			  this.b.translateOf(new Vector_2(-deltaX, -deltaY));
		  
			  this.current=new Point_2(mouseX, mouseY);
		  }
	  }
		
	    /**
	     * Update of the bounding box
	     */    
	    protected void update(double x, double y) {
	    	if (x<xmin)
	    		xmin = x-boundaryThickness;
	    	if (x>xmax)
	    		xmax = x+boundaryThickness;
	    	if (y<ymin)
	    		ymin = y-boundaryThickness;
	    	if (y>ymax)
	    		ymax = y+boundaryThickness;
	    }

	    /**
	     * Update the range of the drawing region (defined by corners points 'a' and 'b')
	     */    
	    public void updateBoundingBox() {
	    	/*
	    	for(Node u: inputGraph.vertices) {
	    		Point_3 p=u.getPoint();
	    		update(p.getX().doubleValue(), p.getY().doubleValue());
	    	}
	    	a=new Point_2(xmin, ymin);
	    	b=new Point_2(xmax, ymax);
	    	*/
	    }
	    
	    /**
	     * Return the current coordinates of the bounding box
	     */    
	    public double[] boundingBox() {
	    	return new double[] {xmin, xmax, ymin, ymax};
	    }

		/**
		 * Return the integer coordinates of a pixel corresponding to a given point
		 * 
		 * Warning: we must take care of the following parameters:
		 * -) the size of the canvas
		 * -) the size of bottom and left panels
		 * -) the negative direction of y-coordinates (in java drawing)
		 */
		public int[] getPoint(Point_2 v) {
			double x=v.getX().doubleValue(); // coordinates of point v
			double y=v.getY().doubleValue();
			double xRange=b.getX().doubleValue()-a.getX().doubleValue(); // width and height of the drawing area
			double yRange=b.getY().doubleValue()-a.getY().doubleValue();
			int i= (int) (this.sizeX*( (x-a.getX().doubleValue()) / xRange )); // scale with respect to the canvas dimension
			int j= (int) (this.sizeY*( (y-a.getY().doubleValue()) / yRange ));
			//i=i+this.horizontalShift;
			j=this.sizeY-j; // y = H - py;
			
			int[] res=new int[]{i, j};
			return res;
		}
		
		  /**
		   * Draw a gray edge (u, v)
		   */
		  public void drawSegment(Point_2 u, Point_2 v) {		  
			int[] min=getPoint(u);
			int[] max=getPoint(v);
		    
			this.stroke(edgeColor, edgeOpacity);
		    this.line(	(float)min[0], (float)min[1], 
		    			(float)max[0], (float)max[1]);
		  }

		  /**
		   * Draw a colored edge (u, v)
		   */
		  public void drawColoredSegment(Point_2 u, Point_2 v, int r, int g, int b) {		  
			int[] min=getPoint(u);
			int[] max=getPoint(v);
		    
			this.stroke(r, g, b, edgeOpacity);
		    this.line(	(float)min[0], (float)min[1], 
		    			(float)max[0], (float)max[1]);
		  }

		  /**
		   * Draw a vertex u on the canvas
		   */
		  public void drawVertex(Node u, double distortion, double maxDistortion) {
			  double ux=u.getPoint().getX().doubleValue();
			  double uy=u.getPoint().getY().doubleValue();
			  double maxValue;
			  
			int[] min=getPoint(new Point_2(ux, uy)); // pixel coordinates of the point in the frame
		    
			//System.out.println("v"+u.index+" dist: "+distortion+" max: "+maxDistortion);
			
			this.stroke(50, 255); // border color
			this.fill(50, 50, 50, 255); // node color
			
			int vertexSize=8; // basic vertex size
			//double growingFactor=1.+(distortion*10.);
			//vertexSize=(int)(3+vertexSize*growingFactor);
			this.ellipse((float)min[0], (float)min[1], vertexSize, vertexSize);
		  }

		  /**
		   * Draw a vertex label on the canvas (close to the node location)
		   */
		  public void drawVertexLabel(Node u) {
			  double ux=u.getPoint().getX().doubleValue();
			  double uy=u.getPoint().getY().doubleValue();			  
			  
			int[] min=getPoint(new Point_2(ux, uy)); // pixel coordinates of the point in the frame
		    			
			String label=this.getVertexLabel(u)+ " weight:" + u.weight; // retrieve the vertex label to show
			
			//this.stroke(edgeColor, edgeOpacity);
			this.fill(200);
			this.rect((float)min[0], (float)min[1], 40, 30); // fill a gray rectangle
			this.fill(0);
			this.text(label, (float)min[0]+5, (float)min[1]+14); // draw the vertex label
		  }

		  /**
		   * Show options on the screen
		   */
		  public void drawOptions() {
			String label="press 'i' or 'o' for zoom\n"; // text to show
			label=label+"press 'c' for one iteration\n";
			label=label+"press 'e' (or 'd') to enable (or to disable) cooling process\n";
			label=label+"press 'l' to output the layout in a 2d frame\n";
			label=label+"use 'left mouse click' to show vertex index\n";
			label=label+"use 'right button button' to drag the layout";
			
			int posX=2;
			int posY=2;
			int textHeight=86;
			
			//this.stroke(edgeColor, edgeOpacity);
			this.fill(240);
			this.rect((float)posX, (float)posY, 380, textHeight); // fill a gray rectangle
			this.fill(0);
			this.text(label, (float)posX+2, (float)posY+10); // draw the text
		  }

		  /**
		   * Select the vertex whose 2d projection is the closest to pixel (i, j)
		   */
		  public Node selectNode(int i, int j) {			  
			  Node result=null;
			  
			  double minDist=40.;
			  for(Node u: inputGraph.vertices) { // iterate over the vertices of g
				  Point_2 p=new Point_2(u.getPoint().getX(), u.getPoint().getY());
				  int[] q=this.getPoint(p);
				  
				  double dist=Math.sqrt((q[0]-i)*(q[0]-i)+(q[1]-j)*(q[1]-j));
				  if(dist<minDist) {
					  minDist=dist;
					  result=u;
				  }
			  }
			  
			  this.selectedNode=result;
			  return result;
		  }
		  
		  /**
		   * Draw the skeleton of a graph in 2D using a Processing frame
		   */
		  public void display2D() {
			  if(this.inputGraph==null)
				  return;
			  AdjacencyListGraph graph=this.inputGraph; // current graph to draw
			  if(graph==null) // if the graph is not defined exit
				  return;
			  
			  this.fill(255,100);
				for(Node u: graph.vertices) { // draw the edges of g
					Point_2 p=new Point_2(u.getPoint().getX(), u.getPoint().getY());
					for(Node v: u.neighbors) {
						if(v!=null && v.index>u.index) { // draw only directed edges (u, v) such that u<v
							Point_2 q=new Point_2(v.getPoint().getX(), v.getPoint().getY());
								this.drawSegment(p, q); // draw a gray edge
						}
					}
				}
				
				for(Node u: graph.vertices) { // finally draw the vertices of g
					this.drawVertex(u, 0, 1.); // color map is not computed
				}

		  }
		  
		  /**
		   * Compute the label of a vertex, from its index, spectral distortion and vertex age
		   */
		  public static String getVertexLabel(Node u) {
		    if (u.getLabel() == null) {
	        String label="v"+u.index;
	        u.setLabel(label);
	        return label;
		    } else {
		      return u.getLabel();}
		  }
		  
			/**
			 * Return an "approximation" (as String) of a given real number (with a given numeric precision)
			 */
			private static String approxNumber(double a, int precision) {
				String format="%."+precision+"f";
				String s=String.format(format,a);
				return s;
			}

}
