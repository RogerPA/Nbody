#include <mpi.h>
#include <stdlib.h>
#include <iostream>
#include <algorithm>
#include <time.h>
#include <iomanip>

using namespace std;


#define DIM 2
				
const double G = 6.673e-11;


int my_id, comm_sz;
MPI_Comm comm;
MPI_Datatype vect_mpi_t;

double** velocities = NULL;


void Fill(double masses[], double** positions, double** local_velocities, int number_particles, int local_number_particles) {
	int part;
	double mass = 5.0e24;
	double gap = 1.0e5;
	double speed = 3.0e4;

	if (my_id == 0) {
		for (part = 0; part < number_particles; part++) {
			masses[part] = mass;
			positions[part][0] = part*gap;
			positions[part][1] = 0.0;
			velocities[part][0] = 0.0;
			if (part % 2 == 0)
				velocities[part][1] = speed;

			else
				velocities[part][1] = -speed;
		}
	}

	MPI_Bcast(masses, number_particles, MPI_DOUBLE, 0, comm);
	MPI_Bcast(positions, number_particles, vect_mpi_t, 0, comm);
	MPI_Scatter(velocities, local_number_particles, vect_mpi_t, local_velocities, local_number_particles, vect_mpi_t, 0, comm);
}

void Print(double masses[], double** positions, double** local_velocities, int number_particles, int local_number_particles) {
	MPI_Gather(local_velocities, local_number_particles, vect_mpi_t, velocities, local_number_particles, vect_mpi_t, 0, comm);

	if (my_id == 0) {
		const int width = 12;
		for (int i = 0; i < number_particles; ++i) {

			cout << i << " " << scientific;
			cout << setw(width) << positions[i][0] << setw(width) << positions[i][1];
			cout << setw(width) << velocities[i][0] << setw(width) << velocities[i][1] << endl;

			//cout << setw(width) << particles[i].force[0] << setw(width) << particles[i].force[1] << endl;

		}
		cout << endl;
	}

}


void ComputeForceBasic(int local_particle, double masses[], double** local_forces, double** position, int number_particles, int local_number_particles) {
	
	int current_particle = my_id * local_number_particles + local_particle;
	double distance_3, x_diff, y_diff, fact;

	local_forces[local_particle][0] = local_forces[local_particle][1] = 0.0;

	for (int i = 0; i < number_particles; ++i) {
		if (i != current_particle) {

			x_diff = position[current_particle][0] - position[i][0];
			y_diff = position[current_particle][1] - position[i][1];
			distance_3 = sqrt(x_diff * x_diff + y_diff * y_diff);
			distance_3 *= distance_3 * distance_3;

			x_diff *= G * masses[current_particle] * masses[i] / distance_3;
			y_diff *= G * masses[current_particle] * masses[i] / distance_3;

			local_forces[local_particle][0] -= x_diff;
			local_forces[local_particle][1] -= y_diff;
		}
	}
}

void UpdateParticles(int local_particle, double masses[], double** local_forces, double** local_position, double** local_velocities, int local_number_particle, double size_step) {

	int current_particle = my_id * local_number_particle + local_particle;
	double step_mass= size_step / masses[current_particle];

	local_position[local_particle][0] += size_step * local_velocities[local_particle][0];
	local_position[local_particle][1] += size_step * local_velocities[local_particle][1];
	local_velocities[local_particle][0] += step_mass * local_forces[local_particle][0];
	local_velocities[local_particle][1] += step_mass * local_forces[local_particle][1];
}

int main(int argc, char* argv[]) {

	int number_particles = 10;
	int local_number_particles;
	int number_step = 4;
	double size_step = 2;
	int current_step;
	int current_particle;


	double* masses;
	double** local_positions;
	double** positions;
	double** local_velocities;
	double** loc_forces;

	MPI_Init(&argc, &argv);
	comm = MPI_COMM_WORLD;
	MPI_Comm_size(comm, &comm_sz);
	MPI_Comm_rank(comm, &my_id);

	local_number_particles = number_particles / comm_sz;
	masses = new double[number_particles];
	positions = new double*[number_particles];
	loc_forces = new double* [local_number_particles];
	local_positions = positions + my_id * local_number_particles;
	local_velocities = new double*[local_number_particles];;

	for (int i = 0; i < number_particles; ++i) {
		positions[i] = new double[DIM];
	}
	for (int i = 0; i < local_number_particles; ++i) {
		loc_forces[i] = new double[DIM];
		local_velocities[i] = new double[DIM];7
	}
	if (my_id == 0) {
		velocities = new double*[number_particles];
		for (int i = 0; i < number_particles; ++i) {
			velocities[i] = new double[DIM];
		}
	}
	MPI_Type_contiguous(DIM, MPI_DOUBLE, &vect_mpi_t);
	MPI_Type_commit(&vect_mpi_t);
	
	
	Fill(masses, positions, local_velocities, number_particles, local_number_particles);

	std::cout.precision(3);

	Print(masses, positions, local_velocities, number_particles, local_number_particles);

	for (current_step = 0; current_step < number_step; ++current_step) {
		for (current_particle = 0; current_particle < local_number_particles; current_particle++) {
			ComputeForceBasic(current_particle, masses, loc_forces, positions, number_particles, local_number_particles);
		}
		for (current_particle = 0; current_particle < local_number_particles; ++current_particle) {
			UpdateParticles(current_particle, masses, loc_forces, local_positions, local_velocities, local_number_particles, size_step);
		}
		MPI_Allgather(MPI_IN_PLACE, local_number_particles, vect_mpi_t, positions, local_number_particles, vect_mpi_t, comm);

		Print(masses, positions, local_velocities, number_particles, local_number_particles);
	}


	MPI_Type_free(&vect_mpi_t);

	MPI_Finalize();

	return EXIT_SUCCESS;
}


