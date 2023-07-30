#include <cmath>
#include <cstdlib>
#include <ctime>
#include<iostream>
#include <vector>
#include <fstream>
#include <sstream>
#include <string>
using namespace std;
struct t_pos
{
	int x;
	int y;
};

double euclidean_distance(double* point1, double* point2, int m_features)
{
	double distance=0;
	for (int i = 0; i < m_features; i++)
	{
		distance += (point1[i] - point2[i])* (point1[i] - point2[i]);
	}
	return std::sqrt(distance);
}

void SOM(t_pos* assignment, double* data, int n_samples,
	int m_features, int height, int width, int max_iter, float lr, float sigma)
{
	std::srand(std::time(nullptr));
	// initialize the weight matrix
	double*** weights = new double** [height];
	for (int i = 0; i < height; i++)
	{
		weights[i] = new double* [width];
		for (int j = 0; j < width; j++)
		{
			weights[i][j] = new double[m_features];
			for (int k = 0; k < m_features; k++)
			{
				weights[i][j][k] = (std::rand() % 100)/1000;
				
			}
		}
		cout << endl;
	}
	for (int t = 0; t < max_iter; t++) {
		// Choose a random sample
		int sample_id = rand() % n_samples;
		double* sample = &data[sample_id * m_features];

		// Find the best matching unit (BMU)
		int bmu_x = 0;
		int bmu_y = 0;
		double min_dist = euclidean_distance(weights[0][0], sample, m_features);
		for (int i = 0; i < height; i++) {
			for (int j = 0; j < width; j++) {
				double dist = euclidean_distance(weights[i][j], sample, m_features);
				if (dist < min_dist) {
					min_dist = dist;
					bmu_x = i;
					bmu_y = j;
				}
			}
		}
		// Update the BMU and its neighbors
		for (int i = 0; i < height; i++) {
			for (int j = 0; j < width; j++) {
				double d = pow(i - bmu_x, 2) + pow(j - bmu_y, 2);
				if (d <= pow(sigma, 2)) {
					for (int k = 0; k < m_features; k++) {
						weights[i][j][k] += lr * (sample[k] - weights[i][j][k]);
					}
				}
			}
		}

		// Decay the learning rate and sigma
		lr *= exp(-t / (double)max_iter);
		sigma *= exp(-t / (double)max_iter);
	}
	
	// Assign each data point to the closest node
	for (int sample_id = 0; sample_id < n_samples; sample_id++) {
		double* sample = &data[sample_id * m_features];

		int bmu_x = 0;
		int bmu_y = 0;
		double min_dist = euclidean_distance(weights[0][0], sample, m_features);
		for (int i = 0; i < height; i++) {
			for (int j = 0; j < width; j++) {
				double dist = euclidean_distance(weights[i][j], sample, m_features);
				if (dist < min_dist) {
					min_dist = dist;
					bmu_x = i;
					bmu_y = j;
				}
			}
		}

		assignment[sample_id].x = bmu_x;
		assignment[sample_id].y = bmu_y;
	}

	// Free memory
	for (int i = 0; i < height; i++) {
		for (int j = 0; j < width; j++) {
			delete[] weights[i][j];
		}
		delete[] weights[i];
	}
	delete[] weights;
}

#define N 180
int main()
{
	std::vector<double> numbers;
	std::ifstream file("BME.txt");  // Open the file

	if (!file)
	{  // Check if the file was opened successfully
		std::cerr << "Unable to open file.txt";
		return 1;   // Call system to stop if can't open file
	}

	double number;
	int count = 0;
	while (file >> number)
	{  // Read numbers from the file
		if (count % 129 != 0)
			numbers.push_back(number);  // Add each number to the vector
		count++;
	}
	file.close();

	std::cout << numbers[0] << std::endl;


		double* data_new = new double[numbers.size()];
		for (int i = 0; i < numbers.size(); i++)
		{
			data_new[i] = numbers[i];
		}
		std::cout << numbers.size()<<" " << data_new[1] << std::endl;
		t_pos assignment[N] ;
		std::cout << "BME:" << std::endl;
		SOM(assignment, data_new, N, 128, 1, 3, 20000, 0.1, 1);
		
		std::cout << "Cluster Result: " << std::endl;
		for (int i = 0; i < N; i++)
			cout << assignment[i].x<< " "<<assignment[i].y<<" |";
		delete[] data_new;

		return 0;
}
