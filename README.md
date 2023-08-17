MOF-Membrane Constructor
The MOF-Membrane Constructor is a Python-based tool for constructing MOF membrane models. This guide will walk you through the installation and execution process.

Installation
Download all four compressed packages in the "MOF_Membrane_Constructor" folder or the scripts compressed in the "source_code.zip" file. Extract them to a directory.

To run the program on Windows (Windows 10/11), click on the Constructor.exe file in the /MOF_Membrane_Constructor folder (extracted from the MOF_Membrane_Constructor.zip.001-004 files). Alternatively, you can run the Constructor.py file from the source code in a Python environment on Windows, Mac OS, or Linux, with the following packages installed: numpy 1.23.0, shapely 1.8.2, and pymatgen 2023.5.31.

To install and run the program using Anaconda, follow these steps:

a) Install shapely: conda install shapely

b) Install pymatgen: conda install pymatgen

c) Run the Constructor: python Constructor.py

Execution
A graphical user interface (GUI) is provided for easy and convenient use of the tool. Follow these steps for a typical construction job:

Click the "Select" button in the "Input Path" panel to choose a directory containing the input structures. Note that only CIF and CAR files are supported as input MOF structures. If using CIF files, we recommend formatting them according to the examples provided in the /Example folder to avoid unwanted output structures.

Define an output directory by clicking the "Select" button in the "Output Path" panel.

Review and adjust the modeling parameters in the middle frame as needed.

Click the "Start" button in the bottom frame to begin the construction process.

Check the "warning.log" file in the output directory for any issues or warnings.

The "Test_for_runtime_and_reliability" directory contains input and output files for evaluating the runtime and reliability of the code. Test results are summarized in the "plane_test.xlsx" and "surface-termination_test.xlsx" files.

For more information on the algorithm and features, please refer to our work (DOI: *****).
