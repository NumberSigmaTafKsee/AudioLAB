% -----------------------------------------------------------------------------
% High-Order Digital Parametric Equalizer Design Toolbox
% -----------------------------------------------------------------------------
% Copyright (c) 2005 by Sophocles J. Orfanidis
% -----------------------------------------------------------------------------
% Address: Sophocles J. Orfanidis                       
%          ECE Department, Rutgers University          
%          94 Brett Road, Piscataway, NJ 08854-8058, USA
%
% Tel:     732-445-5017
% Email:   orfanidi@ece.rutgers.edu
% Date:    June 15, 2005, revised June 19, 2006
% -----------------------------------------------------------------------------
% Reference: Sophocles J. Orfanidis, "High-Order Digital Parametric Equalizer 
%            Design," J. Audio Eng. Soc., vol.53, pp. 1026-1046, November 2005.
% -----------------------------------------------------------------------------
% Web Page: http://www.ece.rutgers.edu/~orfanidi/hpeq
% -----------------------------------------------------------------------------
% tested with MATLAB R11.1 and R14
% -----------------------------------------------------------------------------
%
% Please put all files in the same directory and add it to the matlabpath.
% To get help, just type the function's name at the command prompt.
%
% Equalizer Design Functions
% --------------------------
% hpeq       - high-order digital parametric equalizer design
% hpeq1      - same as hpeq, but uses explicit design equations
% hpeqbw     - bandwidth remapping for high-order digital parametric equalizer
% hpeqex0    - example output from hpeq
% hpeqex1    - example output from hpeq
% hpeqex2    - example output from hpeq
% hpeqord    - order determination of digital parametric equalizer
% bandedge   - calculate left and right bandedge frequencies from bilinear transformation
% blt        - bilinear transformation of analog second-order sections
% octbw      - iterative solution of the octave bandwidth equation
% bandwidth  - calculate bandwidth at any given level 
%
% Elliptic Function Computations
% ------------------------------
% acde       - inverse of cd elliptic function
% asne       - inverse of sn elliptic function
% cde        - cd elliptic function with normalized complex argument
% cne        - cn elliptic function with normalized real argument
% dne        - dn elliptic function with normalized real argument
% sne        - sn elliptic function with normalized complex argument
% ellipdeg   - solves the degree equation (k from N,k1)
% ellipdeg1  - solves the degree equation (k1 from N,k)
% ellipdeg2  - solves the degree equation using nomes
% ellipk     - complete elliptic integral K(k) of first kind
% elliprf    - elliptic rational function
% landen     - Landen transformations of an elliptic modulus
%
% Realizations and Filtering
% --------------------------
% cas2dir    - cascade to direct form coefficients
% dir2decoup - direct form to decoupled form coefficients
% dir2latt   - direct form to normalized lattice form
% dir2state  - direct form to optimum state-space form for second-order filters
% stpeq      - minimum-noise state-space realization of biquadratic digital parametric equalizer
% cascfilt   - filtering in cascade form (uses the built-in function filter)
% decoupfilt - filtering in frequency-shifted 2nd-order cascaded decoupled form
% df2filt    - filtering in frequency-shifted 2nd-order cascaded direct form II
% nlattfilt  - filtering in frequency-shifted 2nd-order cascaded normalized lattice sections
% statefilt  - filtering in frequency-shifted cascaded state-space form
% transpfilt - filtering in frequency-transformed 2nd-order cascaded transposed direct-form-II
%
% Analog Equalizer Design Functions
% ---------------------------------
% bandedge_a - bandedge frequencies for analog equalizer designs
% fresp_a    - frequency response of cascaded analog filter
% hpeq_a     - high-order analog parametric equalizer design
% hpeqbw_a   - remap bandwidth of high-order analog parametric equalizer
% hpeqex1_a  - example output from hpeq_a
% hpeqex2_a  - example output from hpeq_a
% hpeqord_a  - order determination of analog parametric equalizer
%
% Figures from the above Paper
% ----------------------------
% fig10      - Butterworth design example
% fig11      - Chebyshev-1 design example
% fig12      - Chebyshev-2 design example
% fig13      - elliptic design example
% fig14      - designs with common 3-dB widths
% fig15      - time-varying center frequency and bandwidth
% fig16a     - step gain, linear bandwidth
% fig16b     - four-step gain, linear bandwidth
% fig16c     - linear gain, linear bandwidth
%
% Test Scripts
% ------------
% testfilt   - test of seven filtering realizations
% testoctbw  - test of the iterative solution of the bandwidth equation using octbw.m
% teststpeq  - testing of the biquad state-space form
%
% Miscellaneous
% -------------
% abg        - dB to absolute amplitude units
% dbg        - absolute amplitude units to dB
% ellipap2   - substitute for the built-in function ellipap
% ripple     - calculate ripple factors from gains
% xtick      - set xtick marks
% ytick      - set ytick marks
% srem       - symmetrized version of REM (used by ACDE)
% fresp      - frequency response of a cascaded IIR filter at a frequency vector w
% zt         - evaluates z-transform of finite-duration sequence


