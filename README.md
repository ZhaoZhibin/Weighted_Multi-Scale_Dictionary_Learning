# Weighted_Multi-Scale_Dictionary_Learning

[![lic-image]][lic-url]

This repository contains the implementation details of our paper: [Journal of Sound and Vibration]
"[**A weighted multi-scale dictionary learning model and its applications on bearing fault diagnosis**](https://www.sciencedirect.com/science/article/abs/pii/S0022460X19300586)" 
by [Zhibin Zhao](https://zhaozhibin.github.io/). 


## About
Extracting impulsive information under strong background noise and harmonic interference
is a challenging problem for bearing fault diagnosis. Multi-scale transforms have achieved
great success in extracting impulsive feature information, however, how to choose a suitable
transform is a difficult problem, especially in the case of strong noise interference. Therefore,
dictionary learning methods have attracted more and more attention in recent years.
A weighted multi-scale dictionary learning model (WMSDL) is proposed in this paper which
integrates the multi-scale transform and fault information into a unified dictionary learning
model and it successfully overcomes four disadvantages of traditional dictionary learning
algorithms including lacking the multi-scale property; restricting training samples to
local patches; being sensitive to strong harmonic interference; suffering from high computational
complexity. Moreover, algorithmic derivation, computational complexity and parameter
selection are discussed. Finally, The effectiveness of the proposed method is verified by
both the numerical simulations and experiments. Comparisons with other state-of-the-art
methods further demonstrate the superiority of the proposed method.


## Dependencies
- Matlab R2016b
- [[TQWT—toolbox]](http://eeweb.poly.edu/iselesni/TQWT/index.html) from I. W. Selesnick. 
- [[ksvdbox13]](http://www.cs.technion.ac.il/~ronrubin/software.html) from Dr. Ron Rubinstein. 
- [[ompbox10]](http://www.cs.technion.ac.il/~ronrubin/software.html) from Dr. Ron Rubinstein. 

## Pakages

This repository is organized as:
- [funs](https://github.com/ZhaoZhibin/Weighted_Multi-Scale_Dictionary_Learning/tree/master/funs) contains the main functions of the algorithm.
- [util](https://github.com/ZhaoZhibin/Weighted_Multi-Scale_Dictionary_Learning/tree/master/util) contains the extra functions of the test.
- [Results](https://github.com/ZhaoZhibin/Weighted_Multi-Scale_Dictionary_Learning/tree/master/Results) contains the results of the algorithm.
- [tqwt_matlab_toolbox](https://github.com/ZhaoZhibin/Weighted_Multi-Scale_Dictionary_Learning/tree/master/tqwt_matlab_toolbox) contains the TQWT toolbox copied from I. W. Selesnick.
- [ksvdbox13](https://github.com/ZhaoZhibin/Weighted_Multi-Scale_Dictionary_Learning/tree/master/ksvdbox13) contains the KSVD toolbox copied from Dr. Ron Rubinstein.
- [ompbox10](https://github.com/ZhaoZhibin/Weighted_Multi-Scale_Dictionary_Learning/tree/master/ompbox10) contains the OMP toolbox copied from Dr. Ron Rubinstein.
In our implementation, **Matlab R2016b** is used to perform all the experiments.

## Implementation:
Flow the steps presented below:
-  Clone this repository.
```
git clone https://github.com/ZhaoZhibin/Weighted_Multi-Scale_Dictionary_Learning.git
open it with matlab
```
-  Test Simulation: Check the parameters setting of simulation in `Config.m` and run `Test_simulaton.m`. 


## Citation
If you feel our WMSDL is useful for your research, please consider citing our paper: 

```
@article{zhao2019weighted,
  title={A weighted multi-scale dictionary learning model and its applications on bearing fault diagnosis},
  author={Zhao, Zhibin and Qiao, Baijie and Wang, Shibin and Shen, Zhixian and Chen, Xuefeng},
  journal={Journal of Sound and Vibration},
  volume={446},
  pages={429--452},
  year={2019},
  publisher={Elsevier}
}
```
## Contact
- zhibinzhao1993@gmail.com

[lic-image]: https://img.shields.io/aur/license/pac.svg
