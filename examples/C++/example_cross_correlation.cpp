#ifndef SPISEME_HEADER
	#define SPISEME_HEADER
	#include "spiSeMe.h"
#endif
#include <fstream>
#include <sstream>


int main(int argc, char** argv) {

	// Load data
	std::ifstream input_file_A("../data/iei_henon.dat", std::ifstream::in);
	std::ifstream input_file_B("../data/iei_henon_bis.dat", std::ifstream::in);
	if (!input_file_A || !input_file_B) {
		std::cerr << "Error while attempting to open files.\n";
		exit(1);
	}
	double			temp_iei;
	std::string 		line;
	std::istringstream	linestream;
	std::vector<double>	iei_sequence_A;
	std::vector<double>	iei_sequence_B;
	// File readout is rudimentary; it just reads the first column, without any check
	while (std::getline(input_file_A, line)) {
		linestream.str(line);
		linestream >> temp_iei;
		iei_sequence_A.push_back(temp_iei);
		linestream.clear();
	}
	input_file_A.close();
	while (std::getline(input_file_B, line)) {
		linestream.str(line);
		linestream >> temp_iei;
		iei_sequence_B.push_back(temp_iei);
		linestream.clear();
	}
	input_file_B.close();

	// Compute cross correlation of the original sequences
	double bin_width = 0.06, max_lag = 10.0;
	std::vector<double>	C;
	std::vector<double>	lags;
	spiSeMe_return_code 	err;
	err = spiSeMe_event_cross_correlation(C, lags, iei_sequence_A, iei_sequence_B, bin_width, max_lag);
	if (err != SSM_SUCCESS) exit(1);

	// Cross correlation significance estimation
	// --- Surrogate generation
	unsigned int M = 500;
	std::vector< std::vector<double> > iei_surrogates_A;
	std::vector< std::vector<double> > iei_surrogates_B;
	err = spiSeMe_surrogate_jodi(iei_surrogates_A, iei_sequence_A, M);
	if (err != SSM_SUCCESS) exit(1);
	err = spiSeMe_surrogate_jodi(iei_surrogates_B, iei_sequence_B, M);
	if (err != SSM_SUCCESS) exit(1);

	// --- Assessment of cross correlation for each pair of surrogates
	std::vector<double>			C_temp(M, 0.0);
	std::vector< std::vector<double> > 	C_surr;
	for (int i = 0; i < C.size(); i++)
		C_surr.push_back(C_temp);
	for (int i = 0; i < M; i++) {
		err = spiSeMe_event_cross_correlation(C_temp, lags, iei_surrogates_A[i], iei_surrogates_B[i], bin_width, max_lag);
		if (err != SSM_SUCCESS) exit(1);
		for (int j = 0; j < C_surr.size(); j++) {
			C_surr[j][i] = C_temp[j];
		}
	}

	// --- Assessment of significance threshold
	// --- The 1% significance threshold corresponds to the 6th largest value of each bin
	std::vector<double>	C_threshold;
	for (int j = 0; j < C_surr.size(); j++) {
		std::sort(C_surr[j].begin(), C_surr[j].end());
		C_threshold.push_back(C_surr[j][494]);
	}

	// Save cross correlation and significance C_threshold to a file
	std::ofstream output_file("cross_correlation.dat", std::ofstream::out | std::ofstream::trunc);
	output_file << "#lag\tc-corr.\tC_threshold (1% significance level)\n";
	for (int i = 0; i < C_threshold.size(); i++)
		output_file << lags[i] << "\t" << C[i] << "\t" << C_threshold[i] << "\n";

	output_file.close();

	return 0;
}
