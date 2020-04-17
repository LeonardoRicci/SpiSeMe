example: example.cpp
	cd ../../code/C++/; g++ -c spiSeMe_event_autocorrelation.cpp -o ac.o `gsl-config --cflags`
	cd ../../code/C++/; g++ -c spiSeMe_surrogate_jodi.cpp -o jodi.o `gsl-config --cflags`
	cd ../../code/C++/; g++ -c spiSeMe_surrogate_iaaft.cpp -o iaaft.o `gsl-config --cflags`
	cd ../../code/C++/; g++ -c spiSeMe_surrogate_sa.cpp -o sa.o `gsl-config --cflags`
	cd ../../code/C++/; g++ -c spiSeMe_surrogate_dither.cpp -o dither.o `gsl-config --cflags`
	mv ../../code/C++/*.o .
	g++ example.cpp ac.o jodi.o iaaft.o sa.o dither.o `gsl-config --libs` -I../../code/C++/ -o example

example_cross_correlation: example_cross_correlation.cpp
	cd ../../code/C++/; g++ -c spiSeMe_event_cross_correlation.cpp -o xc.o `gsl-config --cflags`
	cd ../../code/C++/; g++ -c spiSeMe_distribution_test.cpp -o kstest.o `gsl-config --cflags`
	cd ../../code/C++/; g++ -c spiSeMe_surrogate_jodi.cpp -o jodi.o `gsl-config --cflags`
	mv ../../code/C++/*.o .
	g++ example_cross_correlation.cpp xc.o jodi.o kstest.o `gsl-config --libs` -I../../code/C++/ -o example_cross_correlation

clean:
	rm -f example example_cross_correlation
	rm -f ac.o xc.o kstest.o jodi.o iaaft.o sa.o dither.o
	rm -f ac_*.dat surr_iei_*.dat cross_correlation.dat
