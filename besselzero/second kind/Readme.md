# Development of Bessel Zeros of the Second Kind
This directory contains the code used to develop the guesses zeros of the Bessel function of the second kind. 

* bruteForceGenerateZeroes.m - Generates zeros using a brutefore search method combined with fzero.
* costFunction.m - helper function used by developGuessEquation.m.
* developGuessEquation.m - Uses the zeros that were saved in zeroes.mat to fit a function that predicts location of the first three zeros given a bessel order.
* zeroes.mat - List of zeros of bessel function of the second kind. 2e6 zeros are contained in the this file corresponding to 200 zeros per order for the first 10,000 orders.

