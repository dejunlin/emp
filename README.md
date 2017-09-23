# emp
Generalized and efficient algorithm for computing multipole energies and gradients based on Cartesian tensors

This is an implementation of an generalized and efficient algorithm to calculate multipolar energies and gradients based on Cartesian tensor. This work is published in

Lin, D. Generalized and efficient algorithm for computing multipole energies and gradients based on Cartesian tensors The Journal of Chemical Physics, 2015, 143, 14115

Please refer to this paper for a thorough discussion on how the algorithm works.

In the current release, an example program is included, where the consistency between numerical gradients and analytic gradients of the Coulomb's kernel 1/r is verified (see the aforementioned paper for details). To build the example program, download and extract the file to a directory and run the following commands:
mkdir build
cd build
ln -sv ../build_release.sh
./build_release.sh
make

Two command line arguments are needed to run the compiled "example" executable, which are the changes in distance (for forces) and quaternion (for torques) for computing the respective numerical gradients. For example,

./example 1e-6 1e-6

This will print out the difference (1st column) as well as the -log10(difference) (2nd column) between the analytic and numerical gradients for each of the pairwise interactions (see the aforementioned paper for details).
