#include <stdlib.h>
#include <iostream>
#include <algorithm>
#include <time.h>
#include <iomanip>
#include <omp.h>


using namespace std;

#define DIM 2

const double G = 6.673e-11;

struct Particle {
	double mass;
	double position[DIM];
	double velocity[DIM];
};

void Fill(Particle* particles, int number_particles) {
	double mass = 5.0e24;
	double gap = 1.0e5;
	double speed = 3.0e4;

	for (int i = 0; i < number_particles; ++i) {
		particles[i].mass = mass;
		particles[i].position[0] = i*gap;
		particles[i].position[1] = 0.0;
		particles[i].velocity[0] = 0.0;
		
		if (i % 2 == 0)
			particles[i].velocity[1] = speed;
		else
			particles[i].velocity[1] = -speed;
	}
}

void Print(Particle* particles, int number_particles) {

	const int width = 12;
	for (int i = 0; i < number_particles; ++i) {

		//cout << i << " " << omp_get_thread_num() <<  scientific;
		cout << i << " " << scientific;
		cout << setw(width) << particles[i].position[0] << setw(width) << particles[i].position[1];
		cout << setw(width) << particles[i].velocity[0] << setw(width) << particles[i].velocity[1] << endl;

		//cout << setw(width) << particles[i].force[0] << setw(width) << particles[i].force[1] << endl;

	}
	cout << endl;
}


void ComputeForceBasic(Particle* particles, double forces[][DIM], int current_particle, int number_particles) {
	forces[current_particle][0] = forces[current_particle][1] = 0.0;
	double distance_3, x_diff, y_diff;
	for (int i = 0; i < number_particles; ++i) {
		if (i != current_particle) {
			x_diff = particles[current_particle].position[0] - particles[i].position[0];
			y_diff = particles[current_particle].position[1] - particles[i].position[1];
			distance_3 = sqrt(x_diff * x_diff + y_diff * y_diff);
			distance_3 *= (distance_3 * distance_3);
			x_diff *= G * particles[i].mass * particles[current_particle].mass / distance_3;
			y_diff *= G * particles[i].mass * particles[current_particle].mass / distance_3;

			forces[current_particle][0] -= x_diff;
			forces[current_particle][1] -= y_diff;
		}
	}

}

void ComputeForceReduced(Particle* particles, double forces[][DIM], int current_particle, int number_particles) {

	double distance_3, x_diff, y_diff;
	for (int i = current_particle + 1; i < number_particles; ++i) {

		x_diff = particles[current_particle].position[0] - particles[i].position[0];
		y_diff = particles[current_particle].position[1] - particles[i].position[1];
		distance_3 = sqrt(x_diff * x_diff + y_diff * y_diff);
		distance_3 *= (distance_3 * distance_3);
		x_diff *= G * particles[i].mass * particles[current_particle].mass / distance_3;
		y_diff *= G * particles[i].mass * particles[current_particle].mass / distance_3;

		forces[current_particle][0] -= x_diff;
		forces[current_particle][1] -= y_diff;
		forces[i][0] += x_diff;
		forces[i][1] += y_diff;
	}

}

void UpdateParticles(Particle* particles, double forces[][DIM], int current_particle, int number_particles, int size_step) {

	double step_mass = size_step / particles[current_particle].mass;

	particles[current_particle].position[0] += size_step * particles[current_particle].velocity[0];
	particles[current_particle].position[1] += size_step * particles[current_particle].velocity[1];
	particles[current_particle].velocity[0] += step_mass * forces[current_particle][0];
	particles[current_particle].velocity[1] += step_mass * forces[current_particle][1];
}

int main(int argc, char* argv[]) {
	
	const int number_threads = 8;
	const int number_particles = 10;
	int number_step = 4;
	double size_step = 2;

	int current_step, current_particle, id_thread;

	Particle particles[number_particles];
	double forces[number_particles][DIM];
	//Reduced
	double local_forces[number_threads * number_particles][DIM];

	std::cout.precision(3);
	srand(time(NULL));

	Fill(particles, number_particles);

	Print(particles, number_particles);

#pragma omp parallel num_threads(number_threads) default(none)	shared(particles, forces, local_forces, size_step, number_particles, number_step) private(current_step, current_particle, id_thread)
	{
		
		id_thread = omp_get_thread_num();
		for (current_step = 0; current_step < number_step; ++current_step) {

//			//Reduced
//#pragma omp for
//			for (current_particle = 0; current_particle < number_threads * number_particles; ++current_particle) {
//				local_forces[current_particle][0] = local_forces[current_particle][1] = 0.0;
//			}
//			//Reduced

//#pragma omp for
#pragma omp for schedule(static, 1)
			// number_step - 1 reduced
			for (current_particle = 0; current_particle < number_particles; ++current_particle) {
				ComputeForceBasic(particles, forces, current_particle, number_particles);
				//ComputeForceReduced(particles, local_forces + id_thread * number_particles, current_particle, number_particles);
			}
//			//Reduced
//# pragma omp for 
//			for (current_particle = 0; current_particle < number_particles; ++current_particle) {
//				forces[current_particle][0] = forces[current_particle][1] = 0.0;
//				for (thread = 0; thread < number_threads; ++thread) {
//					forces[current_particle][0] += local_forces[current_particle + thread * number_particles[0];
//					forces[current_particle][1] += local_forces[current_particle + thread * number_particles][1];
//				}
//			}
//			//Reduced
#pragma omp for
			for (current_particle = 0; current_particle < number_particles; ++current_particle) {
				UpdateParticles(particles, forces, current_particle, number_particles, size_step);
			}
#pragma omp single
			Print(particles, number_particles);
		}
	}

	Print(particles, number_particles);

	return EXIT_SUCCESS;
}