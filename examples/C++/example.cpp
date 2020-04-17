#ifndef SPISEME_HEADER
	#define SPISEME_HEADER
	#include "spiSeMe.h"
#endif
#include <fstream>
#include <sstream>


int main(int argc, char** argv) {

	// Select the sequence in input by checking command line argument
	if (argc < 2) {
		std::cerr << "Please provide a valid label command-line argument...\n";
		std::cerr << "\tExample:\t./example henon\n";
		exit(1);
	}
	std::string label(argv[1]);
	std::string file_name("../data/iei_" + label + ".dat");

	// Load data
	std::ifstream input_file(file_name.c_str(), std::ifstream::in);
	if (!input_file) {
		std::cerr << "Error while attempting to open file " << file_name << ".\n";
		exit(1);
	}
	double			temp_iei;
	std::string 		line;
	std::istringstream	linestream;
	std::vector<double>	iei_sequence;
	while (std::getline(input_file, line)) {	// File readout is rudimentary; it just reads the first column, without any check
		linestream.str(line);
		linestream >> temp_iei;
		iei_sequence.push_back(temp_iei);
		linestream.clear();
	}
	input_file.close();

	// The following parameters are the same as those reported in
	// Chaos 29, 121102 (2019); doi:10.1063/1.5138250
	double bin_width = 0.0, max_lag = 0.0;
	if (!strcmp(label.c_str(), "poisson")) {
		bin_width = 0.03;
		max_lag = 10.5;
	} else if (!strcmp(label.c_str(), "heartbeatlike")) {
		bin_width = 0.03;
		max_lag = 36.0;
	} else if (!strcmp(label.c_str(), "henon")) {
		bin_width = 0.06;
		max_lag = 21.0;
	}

	// Generate surrogate IEIs
	// As an example, three surrogate sequences are generated
	unsigned int M = 3;
	spiSeMe_return_code err;
	std::vector< std::vector<double> > iei_surrogates;
	// --- JODI
	err = spiSeMe_surrogate_jodi(iei_surrogates, iei_sequence, M, true);
	if (err != SSM_SUCCESS) exit(1);
	// --- IAAFT
	err = spiSeMe_surrogate_iaaft(iei_surrogates, iei_sequence, "distribution", M, true);
	if (err != SSM_SUCCESS) exit(1);
	// --- SA
	err = spiSeMe_surrogate_sa(iei_surrogates, iei_sequence, bin_width, max_lag, 0.9, 0.05, 0.64, "max", 200, 50, M, true);
	if (err != SSM_SUCCESS) exit(1);
	// --- DITHERING
	err = spiSeMe_surrogate_dither(iei_surrogates, iei_sequence, "triangular", 0.1, M, true);
	if (err != SSM_SUCCESS) exit(1);

	// Save the surrogate sequences to a file
	std::string output_fname_sequences("surr_iei_" + label + ".dat");
	std::ofstream output_file_sequences(output_fname_sequences, std::ofstream::out | std::ofstream::trunc);
	for (int j = 0; j < iei_surrogates[0].size(); j++) {
		for (int i = 0; i < iei_surrogates.size() - 1; i++) {
			output_file_sequences << iei_surrogates[i][j] << "\t";
		}
		output_file_sequences << iei_surrogates[iei_surrogates.size() - 1][j] << "\n";
	}
	output_file_sequences.close();

	// Compute autocorrelation of original and surrogate sequence
	std::vector<double> lags;
	std::vector<double> autocorr_original;
	std::vector<double> temp_autocorr_surrogate;
	std::vector< std::vector<double> > autocorr_surrogates;
	// --- Autocorrelation of original sequence
	err = spiSeMe_event_autocorrelation(autocorr_original, lags, iei_sequence, bin_width, max_lag);
	if (err != SSM_SUCCESS) exit(1);
	// --- Autocorrelation of surrogate sequences
	for (int i = 0; i < M; i++) {
		err = spiSeMe_event_autocorrelation(temp_autocorr_surrogate, lags, iei_surrogates[i], bin_width, max_lag);
		if (err != SSM_SUCCESS)
			exit(1);

		autocorr_surrogates.push_back(temp_autocorr_surrogate);
	}

	// Save the autocorrelations to a file
	std::string output_fname_autocorrelation("ac_" + label + ".dat");
	std::ofstream output_file_autocorrelation(output_fname_autocorrelation, std::ofstream::out | std::ofstream::trunc);
	output_file_autocorrelation << "#lag\tac_original\tac_surr\n";
	for (int j = 0; j < autocorr_original.size(); j++) {
		output_file_autocorrelation << lags[j] << "\t" << autocorr_original[j];
		for (int i = 0; i < M; i++) {
			output_file_autocorrelation << "\t" << autocorr_surrogates[i][j];
		}
		output_file_autocorrelation << "\n";
	}
	output_file_autocorrelation.close();

	return 0;
}
