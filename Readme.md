
项目源码为张明老师发布于googlecode，众所周知googlecode不方便墙内用户访问，因此在github克隆一份。

可以访问张明老师博客：
http://my.oschina.net/zmjerry/blog

SP++ (Signal Processing in C++) 是一个关于信号处理与数值计算的开源C++程序库，该库提供了信号处理与数值计算中常用算法的C++实现。SP++中所有算法都以C++类模板方法实现，以头文件形式组织而成，所以不需要用户进行本地编译，只要将相关的头文件包含在项目中即可使用。”XXX.h”表示声明文件，”XXX-impl.h”表示对应的实现文件。所有的函数和类均位于名字空间”splab”中，因此使用SP++时要进行命名空间声明：”using namespace splab”。

===============


Signal Processing Library in C++

This is a C++ library for Numerical Computation and Signal Processing. All of the algorithms is implemented by C++ template classes, and organized by ".h" files, so you don't need compile them by yourself.

The SP++ is also published at "Open Source China", the blog publishing URL is: http://my.oschina.net/zmjerry/blog.

The algorithms implemented in SP++ are as follow:

1 Vector Class Template

1.1 Basic Vector Class

1.2 Vector Version for Often Used Functions

1.3 Utilities Functions

1.4 A Simple Timer

2 Matrix Class Template

2.1 Basic Matrix Class

2.2 Matrix Version for Often Used Functions

2.3 Cholesky Decomposition for Real and Complex Matrix

2.4 LU Decomposition for Real and Complex Matrix

2.5 QR Decomposition for Real and Complex Matrix

2.6 SVD Decomposition for Real and Complex Matrix

2.7 Eigenvalue Decomposition for Real and Complex Matrix

2.8 Inversion and Pseudoinversion for Real and Complex Matrix

3 System of Linear Equations

3.1 Common Linear Equations

3.2 Undetermined Linear Equations

3.3 Rank Defect Linear Equations

4 Nonlinear Equation and Equations

4.1 Root of Nonlinear Equation

4.2 Root of Nonlinear Equations

4.3 Romberg Numerical Integration

5 Interpolation and Fitting

5.1 Newton Interpolation

5.2 Cubic Spline Interpolation

5.3 Least Squares Fitting

6 Optimization Method

6.1 Line Searching

6.2 Steepest Descent Method

6.3 Conjugate Gradient Method

6.4 BFGS Method

7 Fourier Transform

7.1 FFT for Signal with Length of 2^n

7.2 FFT for Signal with Arbitrary Length

7.3 A Friendly Used Version of FFT

7.4 C++ Interface for FFTW

7.5 Convolution and Its Fast Algorithm

8 Digital Filter Design

8.1 Widow Functions

8.2 Basic Class for Filter Design

8.3 FIR Digital Filter Design

8.4 IIR Digital Filter Design

9 Random Signal Processing

9.1 Random Number Generator

9.2 Often Used Functions for Probability and Statistics

9.3 Correlation and Its Fast Algorithm

10 Power Spectrum Estimation

10.1 Classical Estimation Methods

10.2 Parameter Estimation Methods

10.3 Eigenanalysis Estimation Methods

11 Adaptive Filters

11.1 Wiener Filter

11.2 Kalman Filter

11.3 LMS Adaptive Filters

11.4 RLS Adaptive Filters

12 Time-Frequency Analysis

12.1 Widow Fourier Transform

12.2 Discrete Gabor Transform

12.3 Wigner-Wille Distribution

13 Wavelet Transform

13.1 Continuous Wavelet Transform

13.2 Dyadic Wavelet Transform

13.3 Discrete Wavelet Transform

14 Searching and Sorting

14.1 Binary Search Tree

14.2 AVL Tree

14.3 Basic Sorting Algorithm

14.4 Huffman Code


