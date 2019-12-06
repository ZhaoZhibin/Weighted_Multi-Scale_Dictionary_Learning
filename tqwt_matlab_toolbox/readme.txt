Tunable Q-Factor Wavelet Transform (TQWT)

Version 1.7

This is a MATLAB implementation of the tunable Q-factor wavelet transform (TQWT) described in the paper.  An implementation in C, which runs faster than the MATLAB code, is available by request.

See the guide: TQWT_guide.pdf for more complete information.

For questions/comments/etc contact selesi@poly.edu

Reference:
'Wavelet Transform with Tunable Q-Factor'
IEEE Trans. on Signal Processing. 59(8):3560-3575, August 2011.
http://eeweb.poly.edu/iselesni/TQWT/
http://taco.poly.edu/selesi/TQWT/
Ivan Selesnick
selesi@poly.edu
Polytechnic Institute of New York University
Brooklyn, NY 11201, USA

--------------------------------------------------

List of functions:

tqwt_radix2.m		Main function for the radix-2 TQWT
itqwt_radix2.m		Inverse radix-2 TQWT
			These functions use radix-2 FFTs only

tqwt.m			TQWT (First form of TQWT described in the TQWT paper)
itqwt.m			Inverse TQWT
			These functions are slower than tqwt_radix2
			The radix-2 versions are preferred for practical use

test_PR_radix2.m	Verify perfect reconstruction properties
test_PR.m

demo1.m			Illustrates the use of functions
demo1_radix2.m
sparsity_demo.m		Reproduces Example 2 in the TQWT paper
resonance_demo.m	Demonstration of resonance decomposition

tqwt_bp.m		Sparse representation using basis pursuit (BP)
tqwt_bpd.m		Sparse approximation using basis pursuit denoising (BPD)
dualQ			Resonance decomposition (BP with dual Q-factors)
dualQd			Resonance decomposition (BPD with dual Q-factors)

PlotFreqResps.m		Plot frequency responses of the TQWT
PlotSubbands.m		Plot subbands computed using TQWT
PlotEnergy		Plot distribution of signal energy over subbands
PlotWavelets.m		Plot wavelets associated with TQWT
ComputeWavelets.m	Compute wavelets associated with TQWT
ComputeNow.m		Compute norms of wavelets
tqwt_fc.m		Center frequencies

afb.m			Analysis filter bank (called by tqwt and tqwt_radix2)
sfb.m			Synthesis filter bank (called by itqwt and itqwt_radix2)
lps.m			Low-pass scaling
next.m			Next power of two
uDFT.m			unitary DFT (normalized DFT)
uDFTinv.m		unitary inverse DFT (normalized inverse DFT)

--------------------------------------------------

Versions

version 1.7		Added mex files for more versions of Matlab

version 1.6		Modify to handle complex data (mfiles)
			Correction to tqwt_bpd mex file
			Miscellaneous other minor changes

version 1.5		April 15, 2011
			Added mex files

version 1.4		April 5, 2011

version 1.3		February 23, 2011

version 1.2		February 22, 2011

version 1.0		January 6, 2011


