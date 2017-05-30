/*
   @author: dephil

   k-dimensional tree structure

   Compilation (command line):
     > javac KDTree.java
     > java KDTree

   Usage (command line):
     > java KDTree

   Note: only executed for testing purposes

 */

// only for plotting
import java.util.*;

import java.awt.Color;
import java.awt.Graphics;
import java.awt.Graphics2D;
import java.awt.Rectangle;
import java.awt.geom.Rectangle2D;
import javax.swing.JFrame;
import javax.swing.JPanel;

public class KDTree extends JPanel {

int K;                    // dimension of the tree
int particlesPerNode;     // maximal number of particles in a leaf node
Node root;                // root node

Particle[] particles;     // particle array


List<double[]> splits;    // (only for plotting)
Particle[] neighbors;     // (only for plotting)
List<double[]> coms;      // (only for plotting)


// only for plotting
final static int XWIN = 1000;
final static int YWIN = 1000;

// Constructor
public KDTree(int particlesPerNode, Particle[] particles) {
        this.K = particles[0].coords().length;      // assuming that all particles have the same dimenions
        this.particlesPerNode = particlesPerNode;
        this.particles = particles;
        initRoot(particles);

        // only for plotting
        this.splits = new ArrayList<double[]>();
        this.coms = new ArrayList<double[]>();
}

// Setter
public void initRoot(Particle[] particles) {
        int numberOfParticles = particles.length;
        int[] min = new int[K];
        int[] max = new int[K];
        for (int d=0; d<K; d++) {
                min[d] = 0;
                max[d] = 0;
                for (int i=1; i<numberOfParticles; i++) {
                        if (particles[min[d]].coords()[d] > particles[i].coords()[d]) {min[d] = i; }
                        if (particles[max[d]].coords()[d] < particles[i].coords()[d]) {max[d] = i; }
                }
        }
        this.root = new Node(0, new double[] {particles[min[0]].coords()[0], particles[min[1]].coords()[1]}, new double[] {particles[max[0]].coords()[0], particles[max[1]].coords()[1]}, 0, numberOfParticles-1);
}

public void spatialBuild(Node base, Particle[] particles) {
        base.branchOffInSpace(particles, base.depth%K);
        splits.add(base.rightChild.rmin);  // only for plotting
        splits.add(base.leftChild.rmax);   // only for plotting
        if (base.leftChild.right-base.leftChild.left >= particlesPerNode)   {
                spatialBuild(base.leftChild, particles);
        }
        if (base.rightChild.right-base.rightChild.left >= particlesPerNode) {
                spatialBuild(base.rightChild, particles);
        }
}

// TODO: particle array is changing constantly which messes up the saved indices when building the tree
public void quickBuild(Node base, Particle[] particles) {
        base.branchOffBalanced(particles, base.depth%K);
        splits.add(base.rightChild.rmin);  // only for plotting
        splits.add(base.leftChild.rmax);   // only for plotting
        if (base.leftChild.right-base.leftChild.left >= particlesPerNode)   {
                //   Partitioner.quickPartition(particles, base.rightChild.left, base.rightChild.right,
                //                              (base.rightChild.right+base.rightChild.left+1)/2, (base.rightChild.depth-1)%K);
                quickBuild(base.leftChild, particles);
        }
        if (base.rightChild.right-base.rightChild.left >= particlesPerNode) {
                // when jumping back up the tree reverse left child ordering before branching off
                if (base.depth>0) {
                        Partitioner.quickPartition(particles, base.leftChild.left, base.leftChild.right,
                                                   (base.leftChild.right+base.leftChild.left+1)/2, (base.leftChild.depth-1)%K);
                }
                quickBuild(base.rightChild, particles);
        }
}

/// Tree methods
public Particle[] findNNs(Node base, Particle p, Particle[] particles, Particle[] neighbors) {
        double key = p.distance2Farthest(neighbors);
        if (base!=null) {
                if (base.distance2Bounds(p.coords())<key) {
                        // if node is close to particle
                        if (base.isLeaf()) {
                                // if node is leaf
                                for (int i=base.left; i<=base.right; i++) {
                                        // if particle i is closer & not the particle in question & not already considered as neighbor
                                        if (particles[i].distance2(p)<key && (p!=particles[i]) && !(Arrays.asList(neighbors).contains(particles[i]))) {
                                                neighbors[(neighbors.length)-1] = particles[i];
                                                key = p.distance2Farthest(neighbors);
                                        }
                                }
                        }
                        else {
                                // if not leaf, go down the tree
                                if (base.leftChild.distance2Bounds(p.coords())<base.rightChild.distance2Bounds(p.coords())) {
                                        neighbors = findNNs(base.leftChild, p, particles, neighbors);
                                        neighbors = findNNs(base.rightChild, p, particles, neighbors);
                                }
                                else {
                                        neighbors = findNNs(base.rightChild, p, particles, neighbors);
                                        neighbors = findNNs(base.leftChild, p, particles, neighbors);
                                }
                        }
                }
        }
        return neighbors;
}

public void leafCOM(Node base, Particle[] particles) {
        // if (base!=null) {
        if (base.isLeaf()) {
                if (base.massInNode(particles)!=0) {
                        base.com = base.calcCOM(particles);
                }
        }
        else {
                leafCOM(base.leftChild, particles);
                leafCOM(base.rightChild, particles);
        }
        // }
}

public void findCOM(Node base, Particle[] particles) {
        if (base!=null) {
                if (base.massInNode(particles)!=0) {
                        base.com = base.calcCOM(particles);
                }
                findCOM(base.leftChild, particles);
                findCOM(base.rightChild, particles);
        }
}

// for plotting - not useful otherwise
public void addCOM(Node base) {
        if (base.isLeaf()) {
                coms.add(base.com);
        }
        else {
                addCOM(base.leftChild);
                addCOM(base.rightChild);
        }
}

// other interesting methods
public int numberOfLeafs(Node base, int nLeafs){
        if (base.isLeaf()) {
                return ++nLeafs;
        } else {
                nLeafs = numberOfLeafs(base.leftChild, nLeafs);
                nLeafs = numberOfLeafs(base.rightChild, nLeafs);
                return nLeafs;
        }
}

public static void printTree(Node base) {
        base.printNode();
        if (base.leftChild!=null) printTree(base.leftChild);
        if (base.rightChild!=null) printTree(base.rightChild);
}

public static void printLeafs(Node base) {
        if (base.isLeaf()) base.printNode();
        if (base.leftChild!=null) printLeafs(base.leftChild);
        if (base.rightChild!=null) printLeafs(base.rightChild);
}


// only for plotting
public void paint(Graphics g) {
        // the particles
        g.setColor(Color.BLACK);
        for (int i=0; i<particles.length; i++) {
                g.fillOval((int)Math.round(XWIN*particles[i].coords()[0]), (int)Math.round(YWIN*particles[i].coords()[1]), 7, 7);
        }
        // draw borders
        g.setColor(Color.GRAY);
        for (int i=0; i<splits.size(); i++) {
                g.drawLine((int)Math.round(splits.get(i)[0]*XWIN), (int)Math.round(splits.get(i)[1]*YWIN), (int)Math.round(splits.get(i+1)[0]*XWIN), (int)Math.round(splits.get(i+1)[1]*YWIN));
                i++;
        }
        // neighbors and target
        g.setColor(Color.RED);
        for (int i=0; i<neighbors.length; i++) {
                g.fillOval((int)Math.round(XWIN*neighbors[i].coords()[0]), (int)Math.round(YWIN*neighbors[i].coords()[1]), 10, 10);
        }
        g.setColor(Color.GREEN);
        g.fillOval((int)Math.round(XWIN*particles[8].coords()[0]), (int)Math.round(YWIN*particles[8].coords()[1]), 10, 10);

        // center of mass
        g.setColor(Color.BLUE);
        for (int i=0; i<coms.size(); i++) {
                if (coms.get(i)!=null) {
                        g.fillOval((int)Math.round(coms.get(i)[0]*XWIN), (int)Math.round(coms.get(i)[1]*YWIN), 4, 4);
                }
        }

        // // particular cell - leftChild
        // if (root.rightChild.leftChild.rightChild.leftChild != null) {
        //         g.setColor(Color.MAGENTA);
        //         for (int i=root.rightChild.leftChild.rightChild.leftChild.left; i<=root.rightChild.leftChild.rightChild.leftChild.right; i++) {
        //                 g.fillOval((int)Math.round(XWIN*particles[i].coords()[0]), (int)Math.round(YWIN*particles[i].coords()[1]), 10, 10);
        //         }
        // }
        //
        // // particular cell - leftChild
        // if (root.leftChild.leftChild.rightChild.leftChild != null) {
        //         g.setColor(Color.CYAN);
        //         for (int i=root.leftChild.leftChild.rightChild.leftChild.left; i<=root.leftChild.leftChild.rightChild.leftChild.right; i++) {
        //                 g.fillOval((int)Math.round(XWIN*particles[i].coords()[0]), (int)Math.round(YWIN*particles[i].coords()[1]), 7, 7);
        //         }
        // }

        // // particular cell - leftChild
        // if (root.rightChild.rightChild.rightChild.leftChild != null) {
        //         g.setColor(Color.ORANGE);
        //         for (int i=root.rightChild.rightChild.rightChild.leftChild.left; i<=root.rightChild.rightChild.rightChild.leftChild.right; i++) {
        //                 g.fillOval((int)Math.round(XWIN*particles[i].coords()[0]), (int)Math.round(YWIN*particles[i].coords()[1]), 5, 5);
        //         }
        // }
        //
        // // particular cell - rightChild
        // if (root.rightChild.rightChild.leftChild.rightChild != null) {
        //         g.setColor(Color.RED);
        //         for (int i=root.rightChild.rightChild.leftChild.rightChild.left; i<=root.rightChild.rightChild.leftChild.rightChild.right; i++) {
        //                 g.fillOval((int)Math.round(XWIN*particles[i].coords()[0]), (int)Math.round(YWIN*particles[i].coords()[1]), 4, 4);
        //         }
        // }
        //
        // // particular cell - rightChild
        // if (root.leftChild.rightChild.rightChild.leftChild != null) {
        //         g.setColor(Color.GREEN);
        //         for (int i=root.leftChild.rightChild.rightChild.leftChild.left; i<=root.leftChild.rightChild.rightChild.leftChild.right; i++) {
        //                 g.fillOval((int)Math.round(XWIN*particles[i].coords()[0]), (int)Math.round(YWIN*particles[i].coords()[1]), 4, 4);
        //         }
        // }
}

/// testing & debugging
public static void main(String[] args) {
        Particle[] particles = Particle.create_2DParticleArray(100);
        KDTree tree = new KDTree(4, particles);
        tree.spatialBuild(tree.root, particles);
        // tree.quickBuild(tree.root, particles);

        tree.neighbors = tree.findNNs(tree.root, particles[8], particles, Particle.create_2DSentinels(4));
        // Particle.print_ParticleArrayCoords(tree.neighbors);

        tree.findCOM(tree.root, particles);
        tree.addCOM(tree.root);


        // for testing
        // KDTree.printTree(tree.root);
        KDTree.printLeafs(tree.root);

        // only for plotting
        JFrame top = new JFrame("KD Tree");
        top.setBounds(10, 10, XWIN, YWIN);
        top.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
        top.getContentPane().add(tree);
        top.setVisible(true);

        tree.paint(top.getGraphics());
}


} /* END KDTREE CLASS ****************************************************** */
