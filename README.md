**Code for 'Unsupervised Co-Learning on $\mathcal{G}$-Manifolds Across Irreducible Representations'** 

**To be appeared in NeurIPS 2019** 

arXiv: https://arxiv.org/abs/1906.02707

## Getting Started

These repository include code examples for the following experiments:

(1) Nearest neighbor search on 2-sphere with G = SO(2);

(2) Nearest neighbor search on 3-sphere with G = SO(3);

(3) Spectral clustering with G = SO(2);

(4) Spectral clustering with G = SO(3);

Code for the cryo-EM image experiment will be updated soon.


## Prerequisites
Code is written in **MATLAB**. To run the codes with 'SO(3)', make sure the _Robotics System Toolbox_ has been installed and authorized.

## Running the Experiments
For example, to reproduce the experiment 'Spectral clustering with G = SO(3)', in Matlab run 
~~~
jobscript_cluster_SO3.m 
~~~
or type 
~~~
jobscript_cluster_SO3
~~~
in the command window.

Each jobscript example has a 'Parameters' section for the hyper-parameter setting. For example:
~~~
n = 1000; % number of data points 
k_max = 6; % maximum frequency
m_k = 10; % eigenvalue truncation for each frequency
id_nn = 50; % number of neighbors for identification
~~~

To adjust the random rewiring probability:
~~~
p_range = [0.16, 0.20, 0.25]; % range of random rewiring probabilities
~~~



## Authors

* **Yifeng Fan** - University of Illinois at Urbana-Champaign

* **Tingran Gao** - University of Chicago (https://gaotingran.com/)

* **Zhizhen Zhao** - University of Illinois at Urbana-Champaign (http://zhizhenz.ece.illinois.edu/)

## Citing
If you find this repository useful for your research, please consider citing our paper:

    @article{fan2019unsupervised,
    title={Unsupervised Co-Learning on $$\backslash$mathcal $\{$G$\}$ $-Manifolds Across Irreducible Representations},
    author={Fan, Yifeng and Gao, Tingran and Zhao, Zhizhen},
    journal={arXiv preprint arXiv:1906.02707},
    year={2019}
    }

## Miscellaneous

Please send any problems/questions you might have to <yifengf2@illinois.edu>.

## License

This project is licensed under the MIT License - see the [LICENSE.md](LICENSE.md) file for details
