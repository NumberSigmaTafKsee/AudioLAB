# Octopus
Octave/Matlab/Scilab

# Matlab Transfer Function
* https://www.mathworks.com/help/control/ug/transfer-functions.html
* https://www.mathworks.com/help/control/ref/tf.html

# Matlab Factor Polynomial
* https://www.mathworks.com/help/symbolic/factor.html

# Matlab Polynomial
* https://www.mathworks.com/help/symbolic/polynomials.html
* https://www.mathworks.com/help/matlab/ref/roots.html


# Octave
* Manual https://docs.octave.org/v7.1.0/index.html#SEC_Contents
* Functions https://octave.sourceforge.io/list_functions.php?sort=alphabetic

# Octave Signal
* Signal Processing https://docs.octave.org/v7.1.0/Signal-Processing.html
* Signal Package https://octave.sourceforge.io/signal/overview.html

# Octave Polynomials
* Findings Roots https://docs.octave.org/v4.2.0/Finding-Roots.html
* Evaluating Polynomial https://docs.octave.org/v4.2.0/Evaluating-Polynomials.html#Evaluating-Polynomials
* Products https://docs.octave.org/v4.2.0/Products-of-Polynomials.html#Products-of-Polynomials
* Calculus https://docs.octave.org/v4.2.0/Derivatives-_002f-Integrals-_002f-Transforms.html#Derivatives-_002f-Integrals-_002f-Transforms
* Interpolation https://docs.octave.org/v4.2.0/Polynomial-Interpolation.html#Polynomial-Interpolation
* Misc https://docs.octave.org/v4.2.0/Miscellaneous-Functions.html#Miscellaneous-Functions


```matlab
syms x
f = 1 + 2*x + 2^x
T = factor(f)

p = sym2poly(f)
zeros(p)
roots(p)
```

# C++ Transfer Function
* https://github.com/borodziejciesla/transfer_function
* https://github.com/KerryL/TransferFunctionPlotter
* https://github.com/DKT19/Symbolic-Transfer-Function
* https://github.com/Tellicious/ArduTF-Library
* https://github.com/MafuraG/CalcTF
* https://github.com/kikuuwe/KTransFunc
* https://github.com/Hammerhead8/Filters
* https://github.com/skrynnyk/Zevaluator
