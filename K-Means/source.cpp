#include<iostream>
#include<vector>
#include<cmath>
#include<algorithm>
#include <fstream>
#include <sstream>
double EuclideanDistance(double* point1, double* point2, int m_features)
{
	double distance=0;
	for (int i = 0; i < m_features; i++)
	{
		distance += (point1[i] - point2[i]) * (point1[i] - point2[i]);
	}
	return std::sqrt(distance);
}

void assign_data_to_centroids(int* assignment, int K, double** centroids, double *data,int n_samples, int m_features)
{
	for (int i = 0; i < n_samples; i++)
	{
		double * current_point = &data[i * m_features];
		int current_cluster = 0;
		double min_distance = EuclideanDistance(current_point, centroids[0], m_features);
		for (int k = 0; k < K; k++)
		{
			double distance = EuclideanDistance(current_point, centroids[k], m_features);
			if (distance <= min_distance)
			{
				current_cluster = k;
				min_distance = distance;
			}
		}
		assignment[i] = current_cluster;
		std::cout << assignment[i];
	}
	std::cout << std::endl;
}

void update_centroids(int * assigment, double** centroids, int K, int n_samples,int m_features, double* data)
{
	int* cluster_size = new int [K] {0};
	for (int i = 0; i < n_samples; i++)
	{
		int cluster = assigment[i];
		cluster_size[cluster]++;
	}
	for (int k = 0; k < K; k++)
	{
		double* centroid = centroids[k];
		std::fill(centroid, centroid + m_features, 0.0);
		for (int i = 0; i < n_samples; i++)
		{
			if (k == assigment[i])
			{
				for (int m = 0; m < m_features; m++)
				{
					centroid[m] += data[i * m_features + m];
				}
			}
		}
		if (cluster_size[k] > 0)
		{
			for (int m = 0; m < m_features; m++)
			{
				centroid[m] /= cluster_size[k];
			}
		}
	}
	delete[] cluster_size;
}

void kmeans(int* assignment, int K, int max_iter, int n_samples, int m_features, double* data)
{
      /*initialize the clusters*/
	double** centroids = new double* [K];
	for (int i = 0; i < K; i++)
	{
		centroids[i] = new double[m_features];
	}
	/* randomly select the initial clusters*/
	std::srand(std::time(nullptr));
	for (int i = 0; i < K; i++)
	{
		int n_cluster = std::rand() % n_samples;
		double* start_point = &data[n_cluster * m_features];
		std::copy(start_point, start_point + m_features, centroids[i]);
	}
	// iterate 
	for (int t = 0; t < max_iter; t++) // iterate for max_iter times
	{
		std::cout << "iterate: " << t << std::endl;
		int* prev_assigment = new int[n_samples];
	    std::copy( assignment, assignment + n_samples, prev_assigment);
	
		assign_data_to_centroids(assignment, K, centroids, data, n_samples, m_features);
		update_centroids(assignment, centroids, K, n_samples, m_features, data);
		bool is_converged = true;
		for (int i = 0; i < n_samples; i++)
		{
			if (assignment[i] != prev_assigment[i])
			{
				is_converged = false;
				break;
			}
		}
		delete[] prev_assigment;
		if (is_converged)
		{
			break;
		}
	}
	for (int i = 0; i < K; i++)
	{
		delete[] centroids[i];
	}
	delete[] centroids;
}
#define N 150
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
		if(count %129 !=0)
		   numbers.push_back(number);  // Add each number to the vector
		count++;
	}
	file.close();




	/*
	std::vector<double> numbers_Iris;
	std::ifstream file_Iris("Iris.txt");  // Open the file
	std::string line;
	if (!file_Iris)
	{  // Check if the file was opened successfully
		std::cerr << "Unable to open file.txt";
		return 1;   // Call system to stop if can't open file
	}

	double number_Iris;
	int count_Iris = 0;
	while (std::getline(file, line)) 
	{
		std::istringstream iss(line);
		char comma;

		for (int i = 0; i < 4; ++i) 
		{
			iss >> number_Iris >> comma;
			numbers_Iris.push_back(number_Iris);
		}
	}

	file_Iris.close();
	std::cout << numbers_Iris[0] << std::endl;
	*/
	std::ifstream file_Iris("iris.txt");  // Open the file
    if (!file_Iris)
	{  // Check if the file was opened successfully
		std::cerr << "Unable to open file.txt";
		return 1;   // Call system to stop if can't open file
	}

	const int NUM_LINES = 150;  // total number of lines in file
	const int NUM_VALUES_PER_LINE = 4;  // number of numeric values per line

	double* data_Iris = new double[NUM_LINES * NUM_VALUES_PER_LINE];
	std::string line;
	int dataIndex = 0;

	while (std::getline(file_Iris, line)) {
		std::istringstream iss(line);
		char comma;

		for (int i = 0; i < NUM_VALUES_PER_LINE; ++i) {
			iss >> data_Iris[dataIndex++] >> comma;
		}
	}
	//std::cout << data_Iris[4] << std::endl;
	file_Iris.close();










	double *data_new = new double[numbers.size()];
	for (int i = 0; i < numbers.size(); i++)
	{
		data_new[i] = numbers[i];
	}
	//std::cout <<numbers.size()<<" " << data_new[1] << std::endl;
	int assignment[N] = { 0 };
	std::cout << "BME:" << std::endl;
	kmeans(assignment, 3, 30, N, 128, data_new);
	std::cout << "Cluster Result: " << std::endl;
	for (int i = 0; i < N; i++)
		std::cout << assignment[i];
	delete[] data_new;
	delete[] data_Iris;
 	return 0;
}