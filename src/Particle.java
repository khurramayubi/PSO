/**
 * Class representing a particle.
 */
public class Particle {

	double[] position; //The position vector of this particle
	double fitness; //The fitness of this particle
	double[] velocity; //The velocity vector of this particle
	double[] personalBest; //Personal best of the particle

	public Particle(double[] position, double[] velocity) {
		this.position = position;
		this.velocity = velocity; 
	}
	
	
}
