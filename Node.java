/*
   @author: dephil

   Nodes for tree structure

   Compilation (command line):
     > javac Node.java
     > java Node
 */

import java.util.Arrays;

public class Node {

int depth;              // depth in the tree
double[] rmin;          // lower spatial bounds
double[] rmax;          // upper spatial bounds
int left;               // lower index bounds in particle array
int right;              // upper index bounds in particle array
double[] com;           // center of mass of node
Node leftChild = null;
Node rightChild = null;
Node parent = null;

/// A bunch of constructors
public Node(int depth) {
        this.depth = depth;
}

public Node(int depth, double[] rmin, double[] rmax) {
        this(depth);
        this.rmin = rmin;
        this.rmax = rmax;
}

public Node(int depth, double[] rmin, double[] rmax, int left, int right) {
        this(depth, rmin, rmax);
        this.left = left;
        this.right = right;
}

public Node(int depth, double[] rmin, double[] rmax, int left, int right, Node parent){
        this(depth, rmin, rmax, left, right);
        this.parent = parent;
}

/// Node branch
public void branchOffInSpace(Particle[] particles, int dimension) {
        double middle = middle(dimension);
        int split = Partitioner.ambiDirPartition(particles, left, right, middle, dimension);
        double[] leftMax = Arrays.copyOf(this.rmax, this.rmin.length);
        leftMax[dimension] = middle;
        double[] rightMin = Arrays.copyOf(this.rmin, this.rmax.length);
        rightMin[dimension] = middle;
        this.leftChild = new Node(this.depth+1, this.rmin, leftMax, left, split-1, this);
        this.rightChild = new Node(this.depth+1, rightMin, this.rmax, split, right, this);
}

public void branchOffBalanced(Particle[] particles, int dimension) {
        int split = (this.right+this.left+1)/2;
        double median = Partitioner.quickPartition(particles, this.left, this.right, split, dimension);
        if ((this.right-this.left+1)%2==0) median = 0.5*(particles[split].coords()[dimension]+Particle.maxCoordInSubarray(particles, this.left, split-1, dimension));
        double[] leftMax = Arrays.copyOf(this.rmax, this.rmax.length);
        leftMax[dimension] = median;
        double[] rightMin = Arrays.copyOf(this.rmin, this.rmin.length);
        rightMin[dimension] = median;
        if ((this.right-this.left+1)%2==1) {
                this.leftChild = new Node(this.depth+1, this.rmin, leftMax, this.left, split, this);
        }
        else {
                this.leftChild = new Node(this.depth+1, this.rmin, leftMax, this.left, split-1, this);
        }
        this.rightChild = new Node(this.depth+1, rightMin, this.rmax, split, this.right, this);
}

/// further Node properties
public boolean isRoot() {
        return (depth == 0) || (parent == null);
}

public boolean isLeaf() {
        return ((leftChild == null) && (rightChild == null));
}

public double size() {
        double size = 0;
        for (int i=0; i<rmax.length; i++) {
                size *= (this.rmax[i] - this.rmin[i]);
        }
        return size;
}

public double middle(int dimension) {
        return 0.5*(this.rmax[dimension]+this.rmin[dimension]);
}

public double[] center(int dimensions) {
        double[] c = new double[dimensions];
        for (int i=0; i<dimensions; i++) {
                c[i] = middle(i);
        }
        return c;
}

public double distance2Bounds(double[] point) { // returns minimal squared distance of a point to the Nodes bounds (and returns 0 if point inside)
        double distance = Double.POSITIVE_INFINITY;
        if (point[1] < this.rmin[1]) {                                  // beneath the node
                if (point[0]<this.rmin[0]) {                              // to the left
                        distance = (this.rmin[0]-point[0])*(this.rmin[0]-point[0]) + (this.rmin[1]-point[1])*(this.rmin[1]-point[1]);
                }
                else if (point[0]>rmax[0]) {                              // to the right
                        distance = (point[0]-this.rmax[0])*(point[0]-this.rmax[0]) + (point[1]-this.rmin[1])*(point[1]-this.rmin[1]);
                }
                else {                                                    // otherwise directly beneath
                        distance = (this.rmin[1]-point[1])*(this.rmin[1]-point[1]);
                }
        }
        else if (point[1] > this.rmax[1]) {                             // above the node
                if (point[0]<this.rmin[0]) {                              // to the left
                        distance = (this.rmin[0]-point[0])*(this.rmin[0]-point[0]) + (this.rmax[1]-point[1])*(this.rmax[1]-point[1]);
                }
                else if (point[0]>rmax[0]) {                              // to the right
                        distance = (point[0]-this.rmax[0])*(point[0]-this.rmax[0]) + (point[1]-this.rmax[1])*(point[1]-this.rmax[1]);
                }
                else {                                                    // otherwise directly above
                        distance = (this.rmax[1]-point[1])*(this.rmax[1]-point[1]);
                }
        }
        else {                                                          // at node level
                if (point[0]<this.rmin[0]) {                              // to the left
                        distance = (this.rmin[0]-point[0])*(this.rmin[0]-point[0]);
                }
                else if (point[0]>this.rmax[0]) {                         // to the right
                        distance = (point[0]-this.rmax[0])*(point[0]-this.rmax[0]);
                }
                else {                                                    // point inside the cell
                        distance = 0;
                }
        }
        return distance;
}


public double massInNode(Particle[] particles) {
        double M = 0;
        if (this.right-this.left>-1) {
                for (int i=this.left; i<=this.right; i++) {
                        M += particles[i].m();
                }
        }
        return M;
}

public double[] calcCOM(Particle[] particles) {
        int dimensions = particles[0].coords().length;  // assuming all particles have the same dimensions
        double[] centerOfMass = new double[dimensions];
        for (int i=0; i<dimensions; i++) {
                centerOfMass[i] = 0;
        }
        double totalMass = massInNode(particles);
        for (int d=0; d<dimensions; d++) {
                for (int i=this.left; i<=this.right; i++) {
                        centerOfMass[d] += particles[i].m()*particles[i].coords()[d];
                }
                centerOfMass[d] /= totalMass;
        }
        return centerOfMass;
}

/// other useful methods
public void printNode() {
        for (int i=0; i<this.depth; i++) {
                System.out.print("-");
        }
        System.out.println(" "+this.left+" "+this.right+"   \t\t-> "+(this.right-this.left+1));
}

} /* END NODE CLASS ******************************************************** */
