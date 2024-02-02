# polylog

 This is Julia code for evaluation of the dilogarithm function. The method is
 based on "The binomial transform of p-recursive sequences and 
 dilogarithm function," [arxiv][def] by Stephanie Harshbarger and Barton Willis. This article was published by _Applications and Applied Mathematics_, in **vol.** 15, **Issue** 2 (December 2020), pp. 1025–1031.

**Features:**

 1. The same code is used for real and complex binary16, binary32, binary64, and bigfloat numbers.

 2. The method uses a series that converges linearly on the entire complex plane. The largest convergence rate is $1/\sqrt{3}$.

 3. To boost accuracy, the method uses Kahan summation. Also, it monitors the accuracy of the summation using a running error bound.



[def]: https://arxiv.org/pdf/1910.06928.pdf