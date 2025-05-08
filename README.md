## The discrete Fourier transform
This repository is a simple implementation of the discrete Fourier transform (DFT) in the C programming language using only standard libraries.

The DFT in this example is run on a sample provided within the code, in this case 18 samples of the function:\
$$f(x) = sin(0.1x/2\pi)+cos(0.3x/2\pi)$$\
\
$$ x_n = {1,2,3... 17,18} \\ N=18 $$

Mathematically, each element k in the DFT is the summation:\
$$F(k)=\sum_{n=0}^{N-1} f(x_n)e^{-i2\pi \frac{k}{N}n $$\
which is implemented by breaking it down into its real and imaginary parts using Euler's identity. The full array of $$F(k)$$ is calculated. The program also outputs $$|F(k)|$$



The inverse of this transform, implemented in the code as iDFT, is as follows:\
$$x_n=\frac{1}{N}\sum_{k=0}^{N-1} F(k)e^{i2\pi \frac{k}{N}n $$\
Where Euler's identity is again used in the implementation to make this calculation purely real. The full sample from $$x_0$$ to $$x_N$$ is reconstructed.
