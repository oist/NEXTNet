#include "stdafx.h"
#include "random.h"
#include "network.h"
#include "analysis.h"

using namespace std;

int program_custered_configuration_model_serrano(int argc, const char * argv[])
{
	std::mt19937 engine;
	const int size = 100;

	// Generate a Poisson graph with the configuration model
	std::poisson_distribution<> poisson(3);
	std::vector<int> degreeList(size,0);
	std::size_t total_degree = 0;
	std::size_t max_degree = 0;
	for (int i = 0; i < size; i++)
	{
		const int k = poisson(engine);
		degreeList[i] = k;
		total_degree += k;
		max_degree = std::max(max_degree, (std::size_t)k);
	}

	// make sure the total degree is even, otherwise no graph can exist
	while (total_degree % 2 == 1) {
		// re-generate a random degree
		const std::size_t i = std::uniform_int_distribution<>(0, size-1)(engine);
		const int d = degreeList[i];
		const int dp = poisson(engine);
		degreeList[i] = dp;
		total_degree += dp - d;
	}

	config_model_clustered_serrano nw(degreeList, 0.5, 0.2, engine);
	std::cerr << "Assortativity: " << assortativity(nw) << std::endl;
	export_adjacency_matrix(nw.adjacencylist, "serrano.txt");
	export_adjacency_dot(nw.adjacencylist, "serrano.dot", false);
	(void)system("dot -Tpdf -Kfdp serrano.dot > serrano.pdf");

	return 0;
}
