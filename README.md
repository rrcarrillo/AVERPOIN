# AVERPOIN (AVErage of squared Residuals for Positive-Input Networks)


## What is it?

AVERPOIN is an application software which evaluates analytically the
representation suitability of a set of neural-network positive inputs.
This suitability refers to the network ability to generate any output
value from any input of the specified set using only positive weights.
Each network input of the set is represented by an array (vector) of
positive values. AVEPOIN considers that the network calculates each output
through the weighted sum of the corresponding input vector values.
Hence AVERPOIN implements an algorithm which receives the set of input
vectors encoded by an m-by-n matrix C. m denotes the number of input
vectors and n is the number of input neurons, thus each row of C
represents an input vector. The output of AVERPOIN is a number between 0
and 1 which encodes the suitability (fitness) of C. A fitness of 0 denotes
the worst possible representation, and 1 a representation for which the
network is able to generate any output value for any input vector.
AVERPOIN is also able to represent graphically the set of inputs so that
its suitability can by visually assessed.
AVERPOIN consists in a set of source code files written in MATLAB
language.
The fitness value returned by AVERPOIN is comparable to the MATLAB
expression rank(C)/size(C,1) but for the case in which the linear
combination of C columns are performed only with positive coefficients.


## Usage

AVERPOIN can be executed by MATLAB and Octave. It has been tested with
MATLAB R2012a and Octave 4.0.0. Considering these versions MATLAB provides
a significantly faster execution and better graphical representation.
The main function for evaluating analytically a cone C is:

> result=int_res(C,plot_flag).

The testing function for evaluating numerically a cone C is:

> result=int_res_num(C,n_pts_dim,plot_flag).

See Documentation section for the function parameter description.

## Documentation

The documentation available as of the date of this release is
included in the .m files. To obtain the usage information of any
function use the "help" command of MATLAB or Octave, especially:

> help int_res

  or

> help int_res_num


## Licensing

Copyright (C) 2017  Richard R. Carrillo (University of Granada)

AVERPOIN is free software: you can redistribute it and/or modify
it under the terms of the GNU Lesser General Public License as
published by the Free Software Foundation, either version 3 of the
License, or (at your option) any later version.

AVERPOIN is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
GNU Lesser General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program (see the files called COPYING.txt and
COPYING_LESSER.txt). If not, see <http://www.gnu.org/licenses/>.


## Contact

* Richard R. Carrillo: rcarrillo .AT. ugr.es
