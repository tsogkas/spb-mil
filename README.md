Basic information
========================

This is the implementation of our symmetry detection system introduced in [1]. It incorporates the Multiple Instance Learning (MIL) framework to train a binary symmetry clasifier on images from the Berkeley Segmentation Dataset (BSDS300) [2]. In the current distribution you will find code for ground-truth construction, training and testing. 
In our work we use the code of the Berkeley pb detector [3], which is publicly available, as well as some parts of the code for the spectral clustering from the Berkeley gpb detector [4]. The image dataset is also publicly available and can be downloaded from the links given below. We also include a publicly available implementation of the skeletonization algorithm in [5] and an implementation of the Lindeberg ridge detector by Iasonas Kokkinos [6].


Compatibility
=========================
The system is implemented in Matlab with some C/C++ functions for efficiency. All mex-files used are already included. Our system works on Linux and Windows, except for the spectral clustering and benchmarking parts, which only work on Linux at the moment. It has been tested on Matlab 2009b - 2011b. There may be some compatibility issues with older versions of Matlab.

For questions concerning the code, feel free to contact Stavros Tsogkas at:
stavros DOT tsogkas AT ecp DOT fr.

References
=========================

[1] Tsogkas, S., Kokkinos I.: Learning-Based Symmetry Detection in Natural Images. ECCV (2012)
 
[2] Martin, D., Fowlkes, C., Tal, D., Malik, J.: A database of human segmented natural
images and its application to evaluating segmentation algorithms and measuring
ecological statistics. In: ICCV (2001).

[3] Martin, D., Fowlkes, C., Malik, J.: Learning to detect natural image boundaries
using local brightness, color, and texture cues. PAMI (2004).

[4] Arbelaez, P., Maire, M., Fowlkes, C., Malik, J.: Contour detection and hierarchical
image segmentation. PAMI (2011)

[5] Telea, A., Van Wijk, J.: An augmented fast marching method for computing skeletons
and centerlines. In: Eurographics (2002).

[6] Kokkinos, I., Maragos, P., Yuille, A.: Bottom-up & top-down object detection using
primal sketch features and graphical models. In: CVPR, vol. 2, pp. 18931900.
IEEE (2006).

Links for Berkeley image dataset and human segmentations:
http://www.eecs.berkeley.edu/Research/Projects/CS/vision/bsds/BSDS300-images.tgz
http://www.eecs.berkeley.edu/Research/Projects/CS/vision/bsds/BSDS300-human.tgz


Usage
==========================
The code for symmetry detection is in the rpb/ directory. You can find separate subdirectories with the code for feature extraction, training, benchmarking etc.

For use on a new example image:
1) Unzip symmetry_1.0 and add it to your working space in Matlab. 
2) The main function for symmetry detection in a new image is rpb.m. Its use is demonstrated in the rpb_demo included. 

For use and training on images from the BSDS300 
1) Add the images and human segmentations in your working path.
2) Edit bsdsRoot.m according to the path of the images and segmentations directory.

Important things to note:
-----------------------------
- In order to train the detector, use the train.m script. You will have to edit the paths of the directories for the ground truth data, the spectral features(if used) and the directory where the resulting weight vector is saved. Specifically you will have to edit paths in train.m and in sampleDetectorMIL.m.

- Spectral feature extraction is time consuming. That is why we prefered to store the extracted features for the desired images and load them when necessary. 

- Training is VERY demanding memory-wise. For training using 100K samples, a computer with at least 16GB ram is necessary. We performed our experiments using 500K samples on servers equipped with 64GB ram. Training time also varies depending on the number of samples, ranging from half an hour to some hours.
