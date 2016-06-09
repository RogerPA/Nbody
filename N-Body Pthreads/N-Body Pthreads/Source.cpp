#define HAVE_STRUCT_TIMESPEC
#include <pthread.h>
#include <stdlib.h>
#include <iostream>
#include <algorithm>
#include <time.h>
#include <iomanip>

using namespace std;

#define DIM 2
const double G = 6.673e-11;

struct Particle {
	double mass;
	double position[DIM];
	double velocity[DIM];
};

const int number_threads = 8;
const int number_particles = 10;
int counter_threads = 0;
int number_step = 4;
double size_step = 2;

int current_step, current_particle;

Particle particles[number_particles];

double forces[number_particles][DIM];
//Reduced
double local_forces[number_threads * number_particles][DIM];

pthread_mutex_t mutex;
pthread_cond_t condition;

void Fill() {
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

void Print() {

	const int width = 12;
	for (int i = 0; i < number_particles; ++i) {
		//    printf("%.3f ", curr[part].m);

		cout << i << " " << scientific;
		cout << setw(width) << particles[i].position[0] << setw(width) << particles[i].position[1];
		cout << setw(width) << particles[i].velocity[0] << setw(width) << particles[i].velocity[1] << endl;

		//cout << setw(width) << particles[i].force[0] << setw(width) << particles[i].force[1] << endl;

	}
	cout << endl;
}


void ComputeForceBasic(int current_particle) {

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

void ComputeForceReduced(int current_particle, double forces[][DIM]) {

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

void UpdateParticles(int current_particle) {

	double step_mass = size_step / particles[current_particle].mass;

	particles[current_particle].position[0] += size_step * particles[current_particle].velocity[0];
	particles[current_particle].position[1] += size_step * particles[current_particle].velocity[1];
	particles[current_particle].velocity[0] += step_mass * forces[current_particle][0];
	particles[current_particle].velocity[1] += step_mass * forces[current_particle][1];
}

void LoopSchedule(int id, int& first, int& last, int& incremental) {

	int quotient = number_particles / number_threads;
	int remainder = number_particles % number_threads;
	int my_iters;
	incremental = 1;
	if (id < remainder) {
		my_iters = quotient + 1;
		first = id * my_iters;
	}
	else {
		my_iters = quotient;
		first = id * my_iters + remainder;
	}
	last = first + my_iters;
}

void WaitCondition() {
	pthread_mutex_lock(&mutex);
	++counter_threads;
	if (counter_threads == number_threads) {
		counter_threads = 0;
		pthread_cond_broadcast(&condition);
	}
	else {
		pthread_cond_wait(&condition, &mutex) != 0;
	}
	pthread_mutex_unlock(&mutex);
}

void* WorkBasic(void* id) {
	int my_id = (int)id;
	int current_step, current_particle;

	for (current_step = 1; current_step <= number_step; ++current_step) {

		for (current_particle = my_id; current_particle < number_particles; current_particle += number_threads) {
			ComputeForceBasic(current_particle);
		}

		WaitCondition();

		for (current_particle = my_id; current_particle < number_particles; current_particle += number_threads) {
			UpdateParticles(current_particle);
		}

		WaitCondition();

		if (my_id == 0) {
			Print();
		}
	}

	return EXIT_SUCCESS;
}

void* WorkReduced(void* id) {
	int my_id = (int)id;
	int current_step, current_particle;
	int first, last, incremental, thread;

	LoopSchedule(my_id, first, last, incremental);

	for (current_step = 1; current_step <= number_step; ++current_step) {

		memset(local_forces + my_id * number_particles, 0, number_particles * sizeof(double[DIM]));

		for (current_particle = first; current_particle < last; current_particle += incremental) {
			ComputeForceReduced(current_particle, local_forces + my_id * number_particles);
		}

		WaitCondition();

		for (current_particle = my_id; current_particle < number_particles; current_particle += number_threads) {
			forces[current_particle][0] = forces[current_particle][1] = 0.0;
			for (thread = 0; thread < number_threads; ++thread) {
				forces[current_particle][0] += local_forces[current_particle + thread * number_particles][0];
				forces[current_particle][1] += local_forces[current_particle + thread * number_particles][1];
			}
		}

		WaitCondition();

		for (current_particle = my_id; current_particle < number_particles; current_particle += number_threads) {
			UpdateParticles(current_particle);
		}

		WaitCondition();

		if (my_id == 0) {
			Print();
		}
	}

	return EXIT_SUCCESS;
}

int main() {

	pthread_t threads[number_threads];

	pthread_mutex_init(&mutex, NULL);
	pthread_cond_init(&condition, NULL);

	std::cout.precision(3);
	srand(time(NULL));

	Fill();

	Print();

	for (int i = 0; i < number_threads; ++i) {
		//pthread_create(&threads[i], NULL, WorkBasic, (void *)i);
		pthread_create(&threads[i], NULL, WorkReduced, (void *)i);
	}

	for (int i = 0; i < number_threads; ++i) {
		pthread_join(threads[i], NULL);
	}

	return EXIT_SUCCESS;
}