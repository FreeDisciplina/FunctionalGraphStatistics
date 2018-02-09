# Experiments on Statistical Properties of the Functional Graph of Random Mappings #

This program perform experimental verification on statistical properties of the functional graph of random mappings.

## Functional Graph of Random Mappings ##

Denote $\mathcal{F}_N$ the set of all mappings from a finite $N$-set into itself. Let $f$ be an arbitrary mapping in $\mathcal{F}_N$. The functional graph of $f$ (denoted by $\mathcal{FG}_{f}$) is a directed graph whose nodes are the elements ${0, \dots, N-1}$ and whose edges are the ordered pairs $\langle x, f(x) \rangle$, for all $x\in \{0, \dots, N-1\}$. It is constructed by successively iterated on $f$. Starting the development by taking an arbitrary node $x_0$ as input and take the output as the input of next iteration on $f$, that is $x_1 = f(x_0), x_2 = f(x_1),\dots,$ before $N$ iterations, a value $x_j$ will equal to one of $x_0, x_1, \dots, x_{j-1}$ ($f$ is a mapping), suppose the one is $x_i$. We say $x_i$ is an `$\alpha$-node` through which the path $x_0\to x_1\to \dots \to x_{i-1}\to  x_{i}$ enters a `cycle` $x_{i}\to x_{i+1}\to \cdots \to x_{j-1} \to x_{i}$. If we consider all possible starting nodes, `paths` exhibit confluence and form into `trees`; trees grafted on cycles form `components`; a collection of components forms a `functional graph` [\[FO90\]](https://link.springer.com/chapter/10.1007/3-540-46885-4_34).

In the functional graph, we call those nodes belongs to a cycle the `cyclic nodes`; those nodes without predecessors (preimages, i.e. $f^{-1}(x) = \emptyset$) the `terminal nodes`; those nodes with at least one predecessor (preimage) the `image nodes`. Seen from an arbitrary node $x_0$, we call the length of the path (measured by the number of edges) starting from $x_0$ and before entering a cycle the `tail length` of $x_0$ and denote by $\lambda(x_0)$; the length of the cycle connected with $x_0$ (measured by the number of edges or nodes) the `cycle length` of $x_0$ and denoted by $\mu(x_0)$; the length of the non repeating trajectory of the node $x_0$ the `rho-length` of $x_0$ and denoted by $\rho(x_0) = \lambda(x_0) + \mu(x_0)$.


The following figure is an illustration of the functional graph of a chopped AES-128 (obtained by fixing an arbitrary key and $122$ bits of the input and take $6$ bits as output, in this case $N = 2^6 = 64$).

<img src="results/Functional_Graph.pdf" />

The structure of the functional graph of random mappings has been studied for a long time.
Lots of parameters have got accurate asymptotic evaluation (refer to [\[FO90\]](https://link.springer.com/chapter/10.1007/3-540-46885-4_34) and [\[FS09\]](http://algo.inria.fr/flajolet/Publications/book.pdf)). These results on statistical properties of the random functional graph can provide various knowledge about the iteration of a random mapping, e.g., the expected number of iterations before one encounter a collision starting from a random node; a quantitative evaluation on the lost entropy of the output when iterating a random function many times. These results have stimulated fruitful results on the cryptanalysis of iterated hash constructions (refer to [\[our Systematization of Knowledge paper\]]()).

## Experimental Verification ##

A practical question is ``For real world pseudo-random mappings designed by cryptographers, how the properties of their functional graph diverse from those statistical properties of the functional graph of random mappings, which is deduced using approaches from analytic combinatorics?'' To answer this question, we performed experiments by simulating a few of $n$-bit random mappings with chopped AES-128 and SM4-128 (obtained by fixing an arbitrary key and $128 - n$ bits of the input and take $n$ bits as the output, $n \in \{12,\ldots, 28\}$). For each $n < 26$, we sample hundreds of random mappings (for $n >= 26$, we can only sample several random mappings due to limited computing resource). We examined the average value, the maximum value, the minimum value and standard deviation respecting to the parameters considered in Theorem 2 in [\[FO90\]](https://link.springer.com/chapter/10.1007/3-540-46885-4_34), e.g., their number of cyclic nodes and their number of $k$-th iterates.

We directly use implementations of AES-128 and SM4-128 included in IPPCP (Cryptography for Intel Intergrated Performance Primitives) in IPP (Intel Integrated Performance Primitives). Thus, one may need to install [IPP](https://software.intel.com/en-us/intel-ipp) and [IPPCP](https://software.intel.com/en-us/get-ipp-cryptography-libraries) to compile and run this experiment.

- This program has been tested on the following platforms:
  -- Windows 10 + VS 2015 + Intel C/C++ compiler 17.0
  -- Red Hat Enterprise Linux Server release 6.9 + g++ (GCC) 5.1.0 + icpc (ICC) 18.0.0

- Compiling:
  -- this is an paralleled implementation, and can be compiled using mpiicpc
  -- four parameters can be defined manually during compiling:
     * `nmin`: the minimum value of $n$ (default 12)
     * `nmax`: the maximum value of $n$ (default 18)
     * `ns`: the number of samples for each $n$-bit random mapping (default 1024)
     * `nt`: the number of mpi processes (default 4)
  -- example:
     * > make nmin=12 nmax=24 ns=768 nt=384

- Running example:
    > mpirun -np 384 ./FG

## Results ##

The following figures show summaries on the results (plotted using [matplotlib](https://matplotlib.org/)):

<img src="results/n12_n25_componentN.pdf" />
<img src="results/n12_n25_cyclicNodeN.pdf" />

<img src="results/n12_n25_tailNodeN.pdf" />
<img src="results/n12_n25_terminalN.pdf" />

<img src="results/n12_n25_imageN.pdf" />
<img src="results/n12_n25_k_thNodeN.pdf" />

## References ##
[FO90] Philippe Flajolet and Andrew M. Odlyzko: Random Mapping Statistics. In Workshop on the Theory and Application of Cryptographic Techniques (EUROCRYPT’89), volume 434 of LNCS, pp. 329–354, Springer, Berlin, Heidelberg, 1990. https://link.springer.com/chapter/10.1007/3-540-46885-4_34

[FS09] Philippe Flajolet and Robert Sedgewick: Analytic Combinatorics. Cambridge University Press, 2009. http://algo.inria.fr/flajolet/Publications/book.pdf
