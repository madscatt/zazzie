#run this file to build the compiled code imported by python under compiledUtils/

cd compiledUtils
python cLoops.py
python cubeScattering.py
cd fortran
python setup_dna_overlap.py build
cd ../
cp fortran/build/*/dna_overlap.so .
cd ../calc_gr/
python setup_gr.py build
cp build/*/gr.so .
