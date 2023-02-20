from setuptools import setup

setup(
	name="rotation_synthesis",
	version="0.0.1",
	install_requires=[
		"numpy",
		"cirq",
		"texttable",
		"pyLIQTR @ git+https://github.com/isi-usc-edu/pyLIQTR.git@main"
		]
	)