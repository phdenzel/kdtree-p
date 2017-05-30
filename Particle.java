/*
   @author: dephil

   Particle data construct

   Compilation (command line):
     > javac Particle.java

   Note: class to be used by Particle N-body code
         can in principle be used for 2D or 3D particles
 */

public class Particle {

public double[] coords;            // particle coords
public double m = 1;               // particle mass
public int numberOfNeighbors = 8;  // number of nearest neighbors
public NearestNeighbors nn;        // nearest neighbors

/// A bunch of constructors
public Particle(double[] point) {
        setCoords(point);
}

public Particle(double[] point, double m) {
        setCoords(point);
        setM(m);
        this.nn = new NearestNeighbors(numberOfNeighbors, point.length);
}

public Particle(double x, double y, double m) { // for 3D particles use constructors above
        this(new double[] {x, y}, m);
}

public Particle(Particle p) {
        setCoords(p.coords());
        setM(p.m());
}

/// Getters
public double[] coords() {
        return this.coords;
}

public double m() {
        return m;
}

public double getX() {
        return this.coords[0];
}

public double getY() {
        return this.coords[1];
}

public double getZ() {
        return this.coords[2];
}

public double distance(Particle otherParticle) { // returns distance to another particle
        int dimension = this.coords().length; // assuming other particle has same dimensions
        double distance = 0;
        double d;
        for (int i=0; i<dimension; i++) {
                d = otherParticle.coords()[i] - this.coords()[i];
                distance += d*d;
        }
        return Math.sqrt(distance);
}

public double distance2(Particle otherParticle) { // returns squared distance to avoid sqrt computations
        int dimension = this.coords().length; // assuming other particle has same dimensions
        double distance = 0;
        double d;
        for (int i=0; i<dimension; i++) {
                d = otherParticle.coords()[i] - this.coords()[i];
                distance += d*d;
        }
        return distance;
}

public double distance(double[] otherPoint) { // returns distance to another point
        int dimension = this.coords().length; // assuming other point has same dimensions
        double distance = 0;
        double d;
        for (int i=0; i<dimension; i++) {
                d = otherPoint[i] - this.coords()[i];
                distance += d*d;
        }
        return Math.sqrt(distance);
}

public double distance2(double[] otherPoint) { // returns squared distance to avoid sqrt computations
        int dimension = this.coords().length; // assuming other point has same dimensions
        double distance = 0;
        double d = 0;
        for (int i=0; i<dimension; i++) {
                d = otherPoint[i] - this.coords()[i];
                distance += d*d;
        }
        return distance;
}

public double distance2Farthest(Particle[] particles) { // returns squared distance to farthest particle in array
        int numberOfParticles = particles.length;
        double dP = 0;
        double dfP = distance2(particles[0]);
        int idfP = 0;
        for (int i=0; i<numberOfParticles; i++) {
                dP = distance2(particles[i]);
                if (dfP < dP) {
                        dfP = dP;
                        idfP = i;
                }
        }
        Partitioner.swap(particles, idfP, particles.length-1);
        return dfP;
}

public int index2Farthest(Particle[] particles) {
        int numberOfParticles = particles.length;
        double dP = 0;
        double dfP = distance2(particles[0]);
        int idfP = 0;
        for (int i=0; i<numberOfParticles; i++) {
                dP = distance2(particles[i]);
                if (dfP < dP) {
                        dfP = dP;
                        idfP = i;
                }
        }
        Partitioner.swap(particles, idfP, particles.length-1);
        return idfP;
}

// more specific distance functions
public double Distance2D(Particle otherParticle) { // 2D distance to another particle
        double dx = otherParticle.getX() - this.getX();
        double dy = otherParticle.getY() - this.getY();
        return Math.sqrt(dx*dx+dy*dy);
}

public double Distance2D(double[] otherPoint) { // 2D distance to another point
        double dx = otherPoint[0] - this.getX();
        double dy = otherPoint[1] - this.getY();
        return Math.sqrt(dx*dx+dy*dy);
}

public double Distance3D(Particle otherParticle) { // 3D distance to another particle
        double dx = otherParticle.getX() - this.getX();
        double dy = otherParticle.getY() - this.getY();
        double dz = otherParticle.getZ() - this.getZ();
        return Math.sqrt(dx*dx+dy*dy+dz*dz);
}

public double Distance3D(double[] otherPoint) { // 3D distance to another point
        double dx = otherPoint[0] - this.getX();
        double dy = otherPoint[1] - this.getY();
        double dz = otherPoint[2] - this.getZ();
        return Math.sqrt(dx*dx+dy*dy+dz*dz);
}

/// Setters
public void setCoords(double[] point) {
        this.coords = point;
}

public void setM(double m) {
        this.m = m;
}

public String coords2String() {
        String particleString = "(";
        for (int i=0; i<coords().length; i++) {
                particleString += Double.toString(coords()[i])+"|";
        }
        particleString = particleString.substring(0, particleString.length() - 2);
        particleString += ")";
        return particleString;
}

/// some particle array methods

public static double minCoordInSubarray(Particle[] particles, int left, int right, int dimension) {
        double min = particles[right].coords()[dimension];
        for (int i=left; i<=right; i++) {
                if (particles[i].coords()[dimension]<min) {
                        min = particles[i].coords()[dimension];
                }
        }
        return min;
}

public static double maxCoordInSubarray(Particle[] particles, int left, int right, int dimension) {
        double max = -1;
        for (int i=left; i<=right; i++) {
                if (particles[i].coords()[dimension]>max) {
                        max = particles[i].coords()[dimension];
                }
        }
        return max;
}

public static void print_ParticleArrayCoords(Particle[] particles){
        for (int i=0; i<particles.length; i++) {
                System.out.print(particles[i].coords2String()+" ");
        }
        System.out.println();
}

public static Particle[] create_2DParticleArray(int numberOfParticles) {
        Particle[] particles = new Particle[numberOfParticles];
        for (int i=0; i<numberOfParticles; i++) {
                particles[i] = new Particle(new double[] {Math.random(), Math.random()});
        }
        return particles;
}

public static Particle[] create_2DSentinels(int numberOfSentinels) {
        Particle[] sentinels = new Particle[numberOfSentinels];
        for (int i=0; i<numberOfSentinels; i++) {
                sentinels[i] = new Particle(new double[] {Double.POSITIVE_INFINITY, Double.POSITIVE_INFINITY});
        }
        return sentinels;
}

public static Particle[] create_3DParticleArray(int numberOfParticles) {
        Particle[] particles = new Particle[numberOfParticles];
        for (int i=0; i<numberOfParticles; i++) {
                particles[i] = new Particle(new double[] {Math.random(), Math.random(), Math.random()});
        }
        return particles;
}

public static Particle[] create_3DSentinels(int numberOfSentinels) {
        Particle[] sentinels = new Particle[numberOfSentinels];
        for (int i=0; i<numberOfSentinels; i++) {
                sentinels[i] = new Particle(new double[] {Double.POSITIVE_INFINITY, Double.POSITIVE_INFINITY, Double.POSITIVE_INFINITY});
        }
        return sentinels;
}

} /* END PARTICLE CLASS ********************************************** */
