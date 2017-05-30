/*
   @author: dephil

   Nearest Neighbors class for tree search

   Compilation (command line):
     > javac NearestNeighbor.java
     > java NearestNeighbor
 */
import java.util.Arrays;


public class NearestNeighbors {

int numberOfNeighbors;  // number of Neighbors
int[] index;            // indices of particle array
double[] distance2;     // squared distances to that neighbor

/// A bunch of constructors
public NearestNeighbors(int numberOfNeighbors, int dimensions) {
        this.numberOfNeighbors = numberOfNeighbors;
        this.index = new int[numberOfNeighbors];
        this.distance2 = new double[numberOfNeighbors];
}

// find nearest neighbors in the kdtree
public Particle[] treeSearch(Node base, Particle p, Particle[] particles, Particle[] neighbors) { // neighbors should first initiated as sentinels
        // NearestNeighbors NN = new NN(neighbors.length, p.coords().length); // assuming all particles have the same dimenions
        double key = p.distance2Farthest(neighbors);
        if (base!=null) {
                if (base.distance2Bounds(p.coords())<key) {
                        // if node is close to particle
                        if (base.isLeaf()) {
                                // if node is leaf
                                for (int i=base.left; i<=base.right; i++) {
                                        if (particles[i].distance2(p)<key && (p!=particles[i]) && !(Arrays.asList(neighbors).contains(particles[i]))) {
                                                // if particle i is closer & not the particle in question & not already considered as neighbor
                                                // TODO: NN.index[neighbors.length-1] = i; NN.distance2[i]
                                                neighbors[(neighbors.length)-1] = particles[i];
                                                key = p.distance2Farthest(neighbors);
                                        }
                                }
                        }
                        else {
                                // if not leaf, go down the tree
                                if (base.leftChild.distance2Bounds(p.coords())<base.rightChild.distance2Bounds(p.coords())) {
                                        neighbors = treeSearch(base.leftChild, p, particles, neighbors);
                                        neighbors = treeSearch(base.rightChild, p, particles, neighbors);
                                }
                                else {
                                        neighbors = treeSearch(base.rightChild, p, particles, neighbors);
                                        neighbors = treeSearch(base.leftChild, p, particles, neighbors);
                                }
                        }
                }
        }
        return neighbors;
}

} /* END NEARESTNEIGHBORS CLASS ******************************************************** */
