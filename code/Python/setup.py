import setuptools

with open("README.md", "r") as fh:
	long_description = fh.read()

setuptools.setup(
	name="spiseme",
	version="1.0.0",
	author="Alessio Perinelli, Michele Castelluzzo, Ludovico Minati and Leonardo Ricci",
	author_email="leonardo.ricci@unitn.it",
	description="Python implementation of four algorithms for surrogate generation of event sequences.",
	long_description=long_description,
	long_description_content_type="text/markdown",
	url="https://github.com/LeonardoRicci/SpiSeMe",
	packages=setuptools.find_packages(),
	classifiers=[
		"Programming Language :: Python :: 3",
		"License :: OSI Approved :: GNU General Public License v3 (GPLv3)",
		"Operating System :: OS Independent",
	],
	python_requires='>=3.5',
	install_requires=[
		'numpy>=1.15'
	]
)
