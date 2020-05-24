%% Main test init
function tests = surrogate_generation_test
	tests = functiontests(localfunctions);
end

%% Test functions
function test_surrogate_jodi_arg(testCase)
	verifyError(testCase, @()spiSeMe_surrogate_jodi(testCase.TestData.data_Z), 'MATLAB:InputParser:ArgumentFailedValidation');
	verifyError(testCase, @()spiSeMe_surrogate_jodi(testCase.TestData.data_A, 'M', 'x'), 'MATLAB:InputParser:ArgumentFailedValidation');
	verifyError(testCase, @()spiSeMe_surrogate_jodi(testCase.TestData.data_A, 'M', 1, 'Verbose', 'x'), 'MATLAB:InputParser:ArgumentFailedValidation');
end
function test_surrogate_jodi_run(testCase)
	y = spiSeMe_surrogate_jodi(testCase.TestData.data_A, 'M', testCase.TestData.M, 'Verbose', false);
	verifyEqual(testCase, size(y,2), testCase.TestData.M);
	verifyEqual(testCase, size(y,1), size(testCase.TestData.data_A,1));
	[p, ~] = spiSeMe_distribution_test(testCase.TestData.data_A, y(:,1));
	verifyEqual(testCase, p, 1.0);
end

function test_surrogate_iaaft_arg(testCase)
	verifyError(testCase, @()spiSeMe_surrogate_iaaft(testCase.TestData.data_Z), 'MATLAB:InputParser:ArgumentFailedValidation');
	verifyError(testCase, @()spiSeMe_surrogate_iaaft(testCase.TestData.data_A, 'exactlyPreserve', 0.1), 'MATLAB:unrecognizedStringChoice');
	verifyError(testCase, @()spiSeMe_surrogate_iaaft(testCase.TestData.data_A, 'exactlyPreserve', 'invalid'), 'MATLAB:unrecognizedStringChoice');
	verifyError(testCase, @()spiSeMe_surrogate_iaaft(testCase.TestData.data_A, 'exactlyPreserve', 'distribution', 'M', 'x'), 'MATLAB:InputParser:ArgumentFailedValidation');
	verifyError(testCase, @()spiSeMe_surrogate_iaaft(testCase.TestData.data_A, 'exactlyPreserve', 'distribution', 'M', 1, 'Verbose', 'x'), 'MATLAB:InputParser:ArgumentFailedValidation');
end
function test_surrogate_iaaft_run(testCase)
	y = spiSeMe_surrogate_iaaft(testCase.TestData.data_A, 'M', testCase.TestData.M, 'Verbose', false);
	verifyEqual(testCase, size(y,2), testCase.TestData.M);
	verifyEqual(testCase, size(y,1), size(testCase.TestData.data_A,1));
	[p, ~] = spiSeMe_distribution_test(testCase.TestData.data_A, y(:,1));
	verifyEqual(testCase, p, 1.0);
end

function test_surrogate_sa_arg(testCase)
	verifyError(testCase, @()spiSeMe_surrogate_sa(testCase.TestData.data_Z, 0.1, 1.0), 'MATLAB:InputParser:ArgumentFailedValidation');
	verifyError(testCase, @()spiSeMe_surrogate_sa(testCase.TestData.data_A, 'x', 1.0), 'MATLAB:InputParser:ArgumentFailedValidation');
	verifyError(testCase, @()spiSeMe_surrogate_sa(testCase.TestData.data_A, 0.1, 'x'), 'MATLAB:InputParser:ArgumentFailedValidation');
 	verifyError(testCase, @()spiSeMe_surrogate_sa(testCase.TestData.data_A, 1.0, 0.1), '');
	verifyError(testCase, @()spiSeMe_surrogate_sa(testCase.TestData.data_A, -0.1, 1.0), '');
	verifyError(testCase, @()spiSeMe_surrogate_sa(testCase.TestData.data_A, 0.1, -1.0), '');
	verifyError(testCase, @()spiSeMe_surrogate_sa(testCase.TestData.data_A, 0.1, 1.0, 'a', 'x'), 'MATLAB:InputParser:ArgumentFailedValidation');
	verifyError(testCase, @()spiSeMe_surrogate_sa(testCase.TestData.data_A, 0.1, 1.0, 'a', 1.1), '');
	verifyError(testCase, @()spiSeMe_surrogate_sa(testCase.TestData.data_A, 0.1, 1.0, 'a', 0.9, 'T', 'x'), 'MATLAB:InputParser:ArgumentFailedValidation');
	verifyError(testCase, @()spiSeMe_surrogate_sa(testCase.TestData.data_A, 0.1, 1.0, 'a', 0.9, 'T', -1.0), '');
	verifyError(testCase, @()spiSeMe_surrogate_sa(testCase.TestData.data_A, 0.1, 1.0, 'a', 0.9, 'T', 0.1, 'C', 1.0, 'costFunction', 'invalid'), 'MATLAB:unrecognizedStringChoice');
	verifyError(testCase, @()spiSeMe_surrogate_sa(testCase.TestData.data_A, 0.1, 1.0, 'a', 0.9, 'T', 0.1, 'C', 1.0, 'costFunction', 'max', 'nTotal', 'x'), 'MATLAB:InputParser:ArgumentFailedValidation');
	verifyError(testCase, @()spiSeMe_surrogate_sa(testCase.TestData.data_A, 0.1, 1.0, 'a', 0.9, 'T', 0.1, 'C', 1.0, 'costFunction', 'max', 'nTotal', 100, 'nSuccessful', 'x'), 'MATLAB:InputParser:ArgumentFailedValidation');
	verifyError(testCase, @()spiSeMe_surrogate_sa(testCase.TestData.data_A, 0.1, 1.0, 'a', 0.9, 'T', 0.1, 'C', 1.0, 'costFunction', 'max', 'nTotal', 100, 'nSuccessful', 50, 'M', 'x'), 'MATLAB:InputParser:ArgumentFailedValidation');
	verifyError(testCase, @()spiSeMe_surrogate_sa(testCase.TestData.data_A, 0.1, 1.0, 'a', 0.9, 'T', 0.1, 'C', 1.0, 'costFunction', 'max', 'nTotal', 100, 'nSuccessful', 50, 'M', 1, 'verbose', 'x'), 'MATLAB:InputParser:ArgumentFailedValidation');
end
function test_surrogate_sa_run(testCase)
	y = spiSeMe_surrogate_sa(testCase.TestData.data_A, 0.1, 1.0, 'a', 0.9, 'T', 0.1, 'C', 1.0, 'costFunction', 'max', 'nTotal', 100, 'nSuccessful', 50, 'M', testCase.TestData.M, 'verbose', false);
	verifyEqual(testCase, size(y,2), testCase.TestData.M);
	verifyEqual(testCase, size(y,1), size(testCase.TestData.data_A,1));
	[p, ~] = spiSeMe_distribution_test(testCase.TestData.data_A, y(:,1));
	verifyEqual(testCase, p, 1.0);
end

function test_surrogate_dither_arg(testCase)
	verifyError(testCase, @()spiSeMe_surrogate_dither(testCase.TestData.data_Z), 'MATLAB:InputParser:ArgumentFailedValidation');
	verifyError(testCase, @()spiSeMe_surrogate_dither(testCase.TestData.data_A, 'ditherDistribution', 3.1), 'MATLAB:unrecognizedStringChoice');
	verifyError(testCase, @()spiSeMe_surrogate_dither(testCase.TestData.data_A, 'ditherDistribution', 'invalid'), 'MATLAB:unrecognizedStringChoice');
	verifyError(testCase, @()spiSeMe_surrogate_dither(testCase.TestData.data_A, 'ditherDistribution', 'uniform', 'distributionParameter', 'x'), 'MATLAB:InputParser:ArgumentFailedValidation');
	verifyError(testCase, @()spiSeMe_surrogate_dither(testCase.TestData.data_A, 'ditherDistribution', 'uniform', 'distributionParameter', 0.1, 'M', 'x'), 'MATLAB:InputParser:ArgumentFailedValidation');
	verifyError(testCase, @()spiSeMe_surrogate_dither(testCase.TestData.data_A, 'ditherDistribution', 'uniform', 'distributionParameter', 0.1, 'M', 1, 'verbose', 'x'), 'MATLAB:InputParser:ArgumentFailedValidation');
end
function test_surrogate_dither_run(testCase)
	y = spiSeMe_surrogate_dither(testCase.TestData.data_A, 'ditherDistribution', 'uniform', 'distributionParameter', 0.1, 'M', testCase.TestData.M, 'verbose', false);
	verifyEqual(testCase, size(y,2), testCase.TestData.M);
	verifyEqual(testCase, size(y,1), size(testCase.TestData.data_A,1));
end

%% Fixtures (setup executed before each function test)
function setup(testCase)
addpath("../");
b = exist('readmatrix');
if (b == 2)
	testCase.TestData.data_A = readmatrix('test_data_A.dat');
else
	testCase.TestData.data_A = dlmread('test_data_A.dat');
end
testCase.TestData.data_Z = [testCase.TestData.data_A, testCase.TestData.data_A];
testCase.TestData.M = 3;
end
