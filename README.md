# MOF-Membrane Constructor
MOF-Membrane Constructor is a Python code that constructs MOF membrane models.
To install and run the code, please follow these steps:


INSTALLATION:

(1) Download all the compressed packages in "program_package" or the scripts in "source_code"; release them in a directory.

(2) Click the Constructor.exe in /dist folder to run this program; This only supports the Windows (Windows 10/11) system. 
Users can also run the Constructor.py of the source code in Python environment to perform this program with the installation of "Shapely" package through the command line in Windows, Mac OS or Linux. 
Specifically, if Anconda software have been installed, users can perform this program in a virtual environments via these steps:
a) conda install shapely ;
b) python Constructor.py ;


EXECUTION:

A GUI is designed for users to quickly and conveniently use this tool. A typical construction job should include the following steps:

(1) Click the "Select" button in the "Input Path" panel to choose a directory that holds the input structures. Note that only CIF and CAR files are supported as the input structures of MOF. Besides, if CIF files are used to describe the input MOFs, we recommend users prepare them in line with the formats of the CIF files in the /Example folder. Otherwise, the program may generate unwanted output structures.

(2) Define a directory to save the output files by clicking the "Select" button in the "Output Path" panel.

(3) Check or set the modeling parameters in the middle frame.

(4) Click the "Start" button in the bottom frame.

(5) Check the "error.log" file in the Output directory for any issues.

For more details about the algorithm and features, please refer to our work (DOI: *****)
