# Inverse Problem of Ultrasound Beamforming with Denoising-Based Regularized Solutions
Sobhan Goudarzi<sup>1</sup>, Adrian Basarab<sup>2</sup>, and Hassan Rivaz<sup>1\
<sup>1</sup> Department of Electrical and Computer Engineering, Concordia University, Montreal, QC, H3G 1M8, Canada.\
<sup>2</sup> Université de Lyon, INSA-Lyon, UCBL, CNRS, Inserm, CREATIS UMR 5220, U1206, Villeurbanne, France.<br><br>
##### Codes used to reproduce the results presented in this [paper](https://arxiv.org/abs/2206.07926), submitted to IEEE Transactions on Ultrasonics, Ferroelectrics, and Frequency Control.
# Abstract
During the past few years, inverse problem formulations of ultrasound beamforming have attracted a growing interest. They usually pose beamforming as a minimization problem of a fidelity term resulting from the measurement model plus a regularization term that enforces a certain class on the resulting image. Herein, we take advantages of alternating direction method of multipliers to propose a flexible framework in which each term is optimized separately. Furthermore, the proposed beamforming formulation is extended to replace the regularization term by a denoising algorithm, based on the recent approaches called plug-and-play (PnP) and regularization by denoising (RED). Such regularizations are shown in this work to better preserve speckle texture, an important feature in ultrasound imaging, than sparsity-based approaches previously proposed in the literature. The efficiency of proposed methods is evaluated on simulations, real phantoms, and *in vivo* data available from a plane-wave imaging challenge in medical ultrasound. Furthermore, a comprehensive comparison with existing ultrasound beamforming methods is also provided. These results show that the RED algorithm gives the best image quality in terms of contrast index while preserving the speckle statistics.
# Requirements
- MATLAB (codes are tested on MATLAB R2021a)
- [PICMUS](https://www.creatis.insa-lyon.fr/Challenge/IEEE_IUS_2016/download) dataset
- [BFGS](https://www.cs.ubc.ca/~schmidtm/Software/minFunc.html) solver
# Usage
1. You need to download the PICMUS dataset and BFGS solver mentioned in the requirements.
2. Once step 1 is completed, the weighting matrix Phi has to be made. To do so, the codes weighting_matrix.m, summation_1.m, and summation_2.m must be run one by one. 
  
# Contact
Sobhan Goudarzi (sobhan.goudarzi@concordia.ca)

# License
[License](https://github.com/Sobhan-Goudarzi/Denoising-Based-Ultrasound-Beamforming/blob/main/LICENSE.txt) for non-commercial use of the software. Please cite the following [paper](https://arxiv.org/abs/2206.07926) when using the codes.
