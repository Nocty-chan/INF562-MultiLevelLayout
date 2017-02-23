import processing.core.PApplet;
import jdg.DrawGraph;
import jdg.graph.AdjacencyListGraph;
import jdg.io.GraphReader;
import jdg.io.GraphReader_MTX;

/**
 * A program for computing network layouts with the "spring embedder" paradigm
 * 
 * This program requires one parameter: the input network, stored in Matrix Market format (.mtx)
 *
 * @author Luca Castelli Aleardi (Ecole Polytechnique, 2017)
 */
public class NetworkLayout extends DrawGraph {
	
	/**
	 * For running the PApplet as Java application
	 */
	public static void main(String args[]) {
		System.out.println("Force-directed layouts (INF562, 2017)");
		if(args.length!=1) {
			System.out.println("Error: wrong arguments, one parameter required: network.mtx");
			//System.out.println("Usage example:\t java -jar NetworkLayout network.mtx");

			System.exit(0);
		}		
		
		String filename=args[0];
		GraphReader reader=new GraphReader_MTX(); // open networks stores in Matrix Market format (.mtx)
		AdjacencyListGraph g=reader.read(filename); // read input network from file
		DrawGraph.inputGraph=g; // set the input network

		PApplet.main(new String[] { "jdg.DrawGraph" }); // start the Processing viewer
	}

}
