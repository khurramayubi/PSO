/**
 * Class representing the PSO Engine. This class implements all the necessary methods for initializing the swarm,
 * updating the velocity and position vectors, determining the fitness of particles and finding the best particle.
 */
public class PSOEngine {

    int numDimensions = 30; //Number of dimensions for problem
    int numParticles = 30; //Number of particles in swarm
    int maxIterations = 10000; //Max number of iterations
    double c1 = 1.496180; //Cognitive coefficient
    double c2 = 1.496180; //Social coefficient
    double w = 0.729844; //Inertia coefficient

    public PSOEngine (int numDimensions, int numParticles, int maxIterations, double c1, double c2, double w ) {
        this.numDimensions = numDimensions;
        this.numParticles = numParticles;
        this.maxIterations = maxIterations;
        this.c1 = c1;
        this.c2 = c2;
        this.w = w;
    }


    /**
     * Method to initialize the particles for PSO
     * @param particles The set of particles to initialize
     */
    public void initParticles(Particle[] particles) {
        //For each particle
        for (int i=0; i<particles.length;i++) {
            double[] positions = new double[numDimensions];
            double[] velocities = new double [numDimensions];
            //For each dimension of the particle assign a random x value [-5.12,5.12] and velocity=0
            for (int j=0; j<numDimensions; j++) {
                positions[j] = ((Math.random()* ((5.12-(-5.12)))) - 5.12);
                velocities[j] = 0;
            }
            //Create the particle
            particles[i] = new Particle(positions, velocities);
            //Set particles personal best to initialized values
            particles[i].personalBest = particles[i].position.clone();
        }
    }

    /**
     * Method to update the velocities vector of a particle
     * @param particle The particle to update the velocity for
     */
    public void updateVelocity(Particle particle, double[] best, double[] r1, double[] r2) {
        //First we clone the velocities, positions, personal and neighbourhood best
        double[] velocities = particle.velocity.clone();
        double[] personalBest = particle.personalBest.clone();
        double[] positions = particle.position.clone();
        double[] bestNeigh = best.clone();

        double[] inertiaTerm = new double[numDimensions];
        double[] difference1 = new double[numDimensions];
        double[] difference2 = new double[numDimensions];

        double[] c1Timesr1 = new double[numDimensions];
        double[] c2Timesr2 = new double[numDimensions];

        double[] cognitiveTerm = new double[numDimensions];
        double[] socialTerm = new double[numDimensions];

        //Calculate inertia component
        for (int i=0; i<numDimensions; i++) {
            inertiaTerm[i] = w*velocities[i];
        }

        //Calculate the cognitive component

        //Calculate personal best - current position
        for (int i=0; i<numDimensions; i++) {
            difference1[i] = personalBest[i] - positions[i];
        }

        //Calculate c1*r1
        for (int i=0; i<numDimensions; i++) {
            c1Timesr1[i] = c1*r1[i];
        }

        //Calculate c1*r1*diff = cognitive term
        for (int i=0; i<numDimensions; i++) {
            cognitiveTerm[i] = c1Timesr1[i]*difference1[i];
        }

        //Calculate the social term

        //Calculate neighbourhood best - current position
        for (int i=0; i<numDimensions; i++) {
            difference2[i] = bestNeigh[i] - positions[i];
        }

        //Calculate c2*r2
        for (int i=0; i<numDimensions; i++) {
            c2Timesr2[i] = c2*r2[i];
        }
        //Calculate c2*r2*diff2 = social component
        for (int i=0; i<numDimensions; i++) {
            socialTerm[i] = c2Timesr2[i]*difference2[i];
        }

        //Update particles velocity at all dimensions
        for (int i=0; i<numDimensions; i++) {
            particle.velocity[i] = inertiaTerm[i]+cognitiveTerm[i]+socialTerm[i];
        }
    }

    /**
     * Method to update the positions vector of a particle
     * @param particle The particle to update the position for
     */

    public void updatePosition(Particle particle) {
        //Since new position is ALWAYS calculated after calculating new velocity, it is okay to just add old position to the current velocity (as velocity would have already been updated).
        for (int i=0; i<numDimensions; i++) {
            particle.position[i] = particle.position[i]+particle.velocity[i];
        }

    }

    /**
     * Method to find the best (fittest) particle from a given set of particles
     * @param particles The collection of particles to determine the best from
     * @return The best (fittest) particle from the collection of particles
     */
    public double[] findBest(Particle[] particles) {
        double[] best = null;
        double bestFitness = Double.MAX_VALUE;
        for(int i=0; i<numParticles; i++) {
            if (evaluateFitness(particles[i].personalBest)<= bestFitness) {
                bestFitness = evaluateFitness(particles[i].personalBest);
                best = particles[i].personalBest;
            }
        }
        return best;
    }

    /**
     * Method to calculate the fitness of a particle using the Rastrigin function
     * @param positions The position vector to evaluate the fitness for
     * @return The fitness of the particle
     */
    public double evaluateFitness(double[] positions) {
        double fitness = 0;
        for (int i=0; i<numDimensions; i++) {
            fitness = fitness + (Math.pow(positions[i],2)-(10*Math.cos(2*Math.PI*positions[i])));
        }

        fitness = fitness + (10*numDimensions);
        return fitness;
    }
}
