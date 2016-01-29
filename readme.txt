To run the three separate parts, run:

python matrixDecomp.py
python iterativeMethod.py
python powerMethod.py

in this directory.

I did use the numpy library. That should be installed in order for this to work properly.

For matrixDecomp.py, up to two files can be passed in. The first is either an Augmented nxn matrix (in which case, that should be the only argument), or just an nxn.
The second argument can be a file containing the 'b' vector.
If no 'b' is found, the zero vector will be used.
If no arguments are provided, then matrixDecomp.py will run for the Pascal matrices (2-12)

example:

python matrixDecomp.py testa.dat testb.dat
python matrixDecomp.py testaug.dat

(Assuming the dat files are in the same directory as partx.py)

For powerMethod, a single file with an n x n matrix can be passed in. To change the initial vector (v), epsilon and the number of iterations, open powerMethod.py in a text editor and change the values at the top of the file. (I will probably reupload with a version that accepts these as arguments...) If a file is not passed in, the program will run the power method on 1000 random matrices.

For matrixDecomp and powerMethod, there can be lot of output... so I recommend running
python partx.py >> out.txt
