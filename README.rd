= DAG kernels for structural RNA analysis -- profile-profile stem kernels --

== Build

You can make the kernel matrix calculator for the profile-profile stem
kernels as follows:

 ./configure &&  make

Several options for the configure script can be specified. You can see
more details with --help option.

The configure script will find required libraries:

* Boost C++ library >= 1.33.1 (http://www.boost.org)
* Vienna RNA packagae >= 1.6  (http://www.tbi.univie.ac.at/~ivo/RNA/)

If you have the MPI system, you may give '--with-mpi' option for the
configure script to build the executable that can run in parallel.


== Learning a Model

Typical usage of the stem kernel is the following command line:

 % stem_kernel_lite -n km.dat +1 pos.fa -1 neg.fa

km.dat is the resulting pre-computed kernel matrix that can be
accepted by LIBSVM (e.g. svm-train -t 4 -b 1 km.dat km.model),
pos.fa and neg.fa are positive sequences and negative sequences,
respectively, written in the FASTA format.
The option -n makes km.dat nomalized to avoid sequence length bias.

The executable can accept the following sequence formats:
* FASTA format,
* CLUSTAL format,
* MAF format.
Furthermore, the filename expansion '*' can be accepted such as:

 % stem_kernel_lite -n km.dat +1 'pos-*.aln' -1 'neg-*.aln'


== Prediction by The Model

 % stem_kernel_lite -n x.dat +1 pos.fa -1 neg.fa --test +1 seq1.fa -1 seq2.fa ....

Before the option --test, you should put the same options as when the
model were learned. After the option --test, you put the class labels
and sequence files to be predicted. Then, you can predict the class of
sequences using LIBSVM (e.g. svm-predict -b 1 x.dat km.model output). 


== References

* Sakakibara, Y., Popendorf, K., Ogawa, N., Asai, K. and Sato, K.: Stem
  kernels for RNA sequence analyses. J Bioinform Comput Biol, 2007, 5,
  1103-1122.
* Sato, K., Kin, T., Asai, K. and Sakakibara, Y.: Directed acyclic graph
  kernels for structural RNA analysis. submitted.

== Contact

((<Kengo SATO|URL:mailto:sato-kengo@aist.go.jp>))
