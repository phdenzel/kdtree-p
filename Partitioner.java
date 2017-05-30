/*
   @author: dephil

   Particle tree partitioner

   Compilation (command line):
     > javac Partitioner.java
     > java Partitioner

   Usage (command line):
     > java Partitioner

   Note: only executed for testing purposes
 */

import java.util.Random;

public class Partitioner {

// one-sided partitioning - returns starting index of second partition
public static int chessboardPartition(Particle[] particles, int left, int right, double value, int dimension) {
    for(int i=left; i<right; i++) {
        if(particles[i].coords()[dimension] < value) swap(particles, i, left++);
    }
    return left;
}

// ambidirectional partitioning - returns starting index of second partition
public static int ambiDirPartition(Particle[] particles, int left, int right, double value, int dimension) {
    while((left<=right) && (particles[left].coords()[dimension] < value)) left++;
    while((left<=right) && (particles[right].coords()[dimension] >= value)) right--;
    while(left<right) {
        swap(particles, left, right);
        left++;
        right--;
        while((left<=right) && (particles[left].coords()[dimension] < value)) left++;
        while((left<=right) && (particles[right].coords()[dimension] >= value)) right--;
    }
    return left;
}

// partitioning method for quickPartition
public static int quickSubStep(Particle[] particles, int left, int right, int pivot, int dimension) {
    double value = particles[pivot].coords()[dimension];
    swap(particles, pivot, right);
    for (int i=left; i<=right; i++) {
        if (particles[i].coords()[dimension]<value) {
            swap(particles, i, left);
            left++;
        }
    }
    swap(particles, right, left);
    return left;
}

// quickSelect paritioning - returns value of second partition's median element
public static double quickPartition(Particle[] particles, int left, int right, int n, int dimension) {
        Random ran = new Random();
        while (right>=left) {
                int pivotIndex = quickSubStep(particles, left, right, ran.nextInt(right-left+1)+left, dimension);
                if (pivotIndex==n) {
                        return particles[pivotIndex].coords()[dimension];
                } else if (pivotIndex<n) {
                        left = pivotIndex + 1;
                } else {
                        right = pivotIndex - 1;
                }
        }
        return particles[right].coords()[dimension];
}

// public static void sortByDistance(Particle particle, Particle[] particles) {
//         return;
// }

public static void swap(Particle[] particles, int i, int j) {
        Particle temp = particles[i];
        particles[i] = particles[j];
        particles[j] = temp;
}

/// only for testing & debugging
public static void print_stats(Particle[] particles, int index, int from, int to, int dimension) {
        for (int i=from; i<to+1; i++) {
                if (i==index) {
                        System.out.println(" ");
                }
                System.out.print(particles[i].coords()[dimension]+" ");
        }
        System.out.println(" ");
        System.out.println("Partition index: "+index);
}

public static void main(String[] args) {
        int noP, dim, index=-1;
        if(args.length==0) {
                noP = 8;
                dim = 0;
        } else {
                noP = Integer.parseInt(args[0]);
                if (args[1].equals("x") || args[1].equals("X") || args[1].equals("0")) {
                        dim = 0;
                }
                else {
                        dim = 1;
                }
        }
        Particle[] particles = Particle.create_2DParticleArray(noP);
        // index = Partitioner.ambiDirPartition(particles, 0, noP-1, 0.5, dim);
        Partitioner.print_stats(particles, 0, 0, noP-1, dim);
        double median = Partitioner.quickPartition(particles, 0, noP-1, noP/2, dim);
        if (index==-1) Partitioner.print_stats(particles, 0, 0, noP-1, dim);
        else Partitioner.print_stats(particles, index, 0, noP-1, dim);
        System.out.println(Particle.maxCoordInSubarray(particles, 0, noP/2-1, dim));
        System.out.println(median);
}

} /* END PARTITIONER CLASS ************************************************** */
