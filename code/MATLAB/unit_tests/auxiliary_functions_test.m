%% Main test init
function tests = auxiliary_functions_test
	tests = functiontests(localfunctions);
end

%% Test functions
function test_aux_distrtest_arg(testCase)
	verifyError(testCase, @()spiSeMe_distribution_test(testCase.TestData.data_Z, testCase.TestData.data_A), 'MATLAB:InputParser:ArgumentFailedValidation');
	verifyError(testCase, @()spiSeMe_distribution_test(testCase.TestData.data_A, testCase.TestData.data_Z), 'MATLAB:InputParser:ArgumentFailedValidation');
end
function test_aux_distrtest_run(testCase)
	[p, d] = spiSeMe_distribution_test(testCase.TestData.data_A, testCase.TestData.data_A);
	verifyEqual(testCase, p, 1.0);
	verifyEqual(testCase, d, 0.0);
	[p, d] = spiSeMe_distribution_test(testCase.TestData.data_A, testCase.TestData.data_B);
	verifyEqual(testCase, p, 0.5, 'AbsTol', 0.1);
	verifyEqual(testCase, d, 0.04, 'AbsTol', 0.01);
end

function test_aux_autocorr_arg(testCase)
	verifyError(testCase, @()spiSeMe_event_autocorrelation(testCase.TestData.data_Z, 0.1, 1.0), 'MATLAB:InputParser:ArgumentFailedValidation');
	verifyError(testCase, @()spiSeMe_event_autocorrelation(testCase.TestData.data_A, 'x', 1.0), 'MATLAB:InputParser:ArgumentFailedValidation');
	verifyError(testCase, @()spiSeMe_event_autocorrelation(testCase.TestData.data_A, 0.1, 'x'), 'MATLAB:InputParser:ArgumentFailedValidation');
	verifyError(testCase, @()spiSeMe_event_autocorrelation(testCase.TestData.data_A, -0.1, 1.0), '');
	verifyError(testCase, @()spiSeMe_event_autocorrelation(testCase.TestData.data_A, 0.1, -1.0), '');
	verifyError(testCase, @()spiSeMe_event_autocorrelation(testCase.TestData.data_A, 1.0, 0.1), '');
end
function test_aux_autocorr_run(testCase)
	[A, lags] = spiSeMe_event_autocorrelation(testCase.TestData.data_A, 0.1, 1.0);
	verifyEqual(testCase, isvector(A), true);
	verifyEqual(testCase, isvector(lags), true);
	verifyEqual(testCase, size(A,1), 10);
	verifyEqual(testCase, size(lags,1), size(A,1));
end

function test_aux_crosscorr_arg(testCase)
	verifyError(testCase, @()spiSeMe_event_cross_correlation(testCase.TestData.data_Z, testCase.TestData.data_A, 0.1, 1.0), 'MATLAB:InputParser:ArgumentFailedValidation');
	verifyError(testCase, @()spiSeMe_event_cross_correlation(testCase.TestData.data_A, testCase.TestData.data_Z, 0.1, 1.0), 'MATLAB:InputParser:ArgumentFailedValidation');
	verifyError(testCase, @()spiSeMe_event_cross_correlation(testCase.TestData.data_A, testCase.TestData.data_B, 'x', 1.0), 'MATLAB:InputParser:ArgumentFailedValidation');
	verifyError(testCase, @()spiSeMe_event_cross_correlation(testCase.TestData.data_A, testCase.TestData.data_B, 0.1, 'x'), 'MATLAB:InputParser:ArgumentFailedValidation');
	verifyError(testCase, @()spiSeMe_event_cross_correlation(testCase.TestData.data_A, testCase.TestData.data_B, -0.1, 1.0), '');
	verifyError(testCase, @()spiSeMe_event_cross_correlation(testCase.TestData.data_A, testCase.TestData.data_B, 0.1, -1.0), '');
	verifyError(testCase, @()spiSeMe_event_cross_correlation(testCase.TestData.data_A, testCase.TestData.data_B, 1.0, 0.1), '');
end
function test_aux_crosscorr_run(testCase)
	[C, lags] = spiSeMe_event_cross_correlation(testCase.TestData.data_A, testCase.TestData.data_B, 0.1, 1.0);
	verifyEqual(testCase, isvector(C), true);
	verifyEqual(testCase, isvector(lags), true);
	verifyEqual(testCase, size(C,1), 21);
	verifyEqual(testCase, size(lags,1), size(C,1));
end

%% Fixtures (setup executed before each function test)
function setup(testCase)
addpath("../");
b = exist('readmatrix');
if (b == 2)
	testCase.TestData.data_A = readmatrix('test_data_A.dat');
	testCase.TestData.data_B = readmatrix('test_data_B.dat');
else
	testCase.TestData.data_A = dlmread('test_data_A.dat');
	testCase.TestData.data_B = dlmread('test_data_B.dat');
end
testCase.TestData.data_Z = [testCase.TestData.data_A, testCase.TestData.data_A];
testCase.TestData.M = 3;
end
