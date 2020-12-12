/**
 * Class implementing the PSO algorithm. 
 */
public class PSOimplementation {

	public final int numDimensions = 30; //Number of dimensions for problem
	public final int numParticles = 30; //Number of particles in swarm
	public final int maxIterations = 10000; //Max number of iterations
	public final double c1 = 1.496180; //Cognitive coefficient
	public final double c2 = 1.496180; //Social coefficient
	public final double w = 0.729844; //Inertia coefficient
	public  double[] r1; //Random vector 1
	public  double[] r2;  //Random vector 2
	public double[] best;
	Particle[] particles; //Array to hold all particles
	
	public PSOimplementation() {
		//PSO algorithm

		particles = new Particle[numParticles];
		PSOEngine PSO = new PSOEngine(numDimensions, numParticles, maxIterations, c1, c2, w);

		//Initialize particles
		PSO.initParticles(particles);

		//PSO loop
		int numIter = 0;
		while (numIter<maxIterations) {
			// Evaluate fitness of each particle
			for (int i=0; i<numParticles; i++) {
				particles[i].fitness = PSO.evaluateFitness(particles[i].position);

				//update personal best position 
				if (particles[i].fitness <= PSO.evaluateFitness(particles[i].personalBest)) {
					particles[i].personalBest = particles[i].position.clone();
				}
			}
			//Find best particle in set
			best = PSO.findBest(particles);

			//Initialize the random vectors for updates
			r1 = new double[numDimensions];
			r2 = new double[numDimensions];
			for (int i=0; i<numDimensions; i++) {
				r1[i] = Math.random();
				r2[i] = Math.random();
			}

			//Update the velocity and position vectors
			for (int i=0; i<numParticles;i++) {
				PSO.updateVelocity(particles[i], best, r1, r2);
				PSO.updatePosition(particles[i]);
			}
			numIter++;
		}
		//Print the best solution
		print((best));
		System.out.println(PSO.evaluateFitness(best));
	}


	/**
	 * Helped method to print an array as a vector
	 * @param a The given 1-D array
	 */
	public void print (double[] a) {
		System.out.print("< ");
		for (int i=0; i<a.length; i++) {
			System.out.print(a[i]  + " ");
		}
		System.out.println(" > ");

	}
	
	public static void main(String[] args) {
		PSOimplementation p = new PSOimplementation(); 
	}

}
