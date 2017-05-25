# spb-mil
Use Multiple Instance Learning (MIL) to train a binary classifier of local symmetries on images from the Berkeley Segmentation Dataset (BSDS300). 

[Learning-Based Symmetry Detection for Natural Images](http://tsogkas.github.io/publications/symmetry-mil/tsogkas2012learning.pdf)  
[Stavros Tsogkas](http://tsogkas.github.io/), [Iasonas Kokkinos](http://www0.cs.ucl.ac.uk/staff/I.Kokkinos/index.html)  
In ECCV, 2012.


## License

This code is released under the MIT License (refer to the LICENSE file for details).

## Contents
1. [Requirements: software](#requirements-software)
2. [Requirements: hardware](#requirements-hardware)
3. [Directory structure](#directory-structure)
4. [Setup](#setup)
5. [Using the code](#using-the-code)
6. [Citation](#citation)
7. [References](#references)


## Requirements: software
* Linux OS (tested on 16.04).
* A recent version of MATLAB. All our experiments were performed using MATLAB R2016a but previous versions of our
code have been successfully tested on MATLAB R2009b-R2011b.
* [`matlab-utils`](https://github.com/tsogkas/matlab-utils).

## Requirements: hardware

You can run our code in any modern computer (desktop or laptop). `spbMIL` runs at ~15sec for a 481x321 image on a modern desktop CPU.

## Directory structure
Generally:
* All data should go under `data/`.
* All external code should go under `external/`.
* `spb-mil` results, models etc should go under `output/`.
* Specific results and plots are saved in the respective directories, e.g.:
  - trained models and medial point detection results are saved in `output/models/`. 
  - medial points detection plots are saved in `output/plots/`.

Feel free to change the paths in `setPaths.m` and use symbolic links to change directory hierarchy to your preference.

## Setup

1. Clone the `spb-mil` repository: `git clone git@github.com:tsogkas/spb-mil.git` and add `cpp`, `external` and `util` to your working path.
2. Clone the `matlab-utils` repository: `git clone git@github.com:tsogkas/matlab-utils.git` and add it to your working path.
3. Create folders `output/`, `data/`.
4. Download the BSDS300 [images](http://www.eecs.berkeley.edu/Research/Projects/CS/vision/bsds/BSDS300-images.tgz) and [human annotations](http://www.eecs.berkeley.edu/Research/Projects/CS/vision/bsds/BSDS300-human.tgz)
and add them in `data/BSDS300`. If you want to use the newer version of the dataset, you can download the [BSDS500](http://www.eecs.berkeley.edu/Research/Projects/CS/vision/grouping/BSR/BSR_bsds500.tgz) dataset and benchmark code and extract it in `data/`. 
5. Download the [SYMMAX300](http://tsogkas.github.io/publications/symmetry-mil/SYMMAX300.zip) dataset.

*NOTE: Please edit the `setPaths` function accordingly so that the paths reflect your directory structure.*

## Using the code

### Training
You can train our MIL-detector using the `trainMIL` function. For example, to train a model on the _train_ set of BSDS, using color and texture features, with 1000 training samples per image, use: 

	trainMIL('trainSet','train','featureSet','color','nSamplesPerImage',1000);

### Testing
You can run performance evaluation tests on the _val_ subset using the command:

	model = testSPB('modelPath', 'dataset','SYMMAX300', 'testSet','val');

Performance statistics are contained in the `model.SYMMAX300.val.stats` struct.

### Running the detector on an input image

	img = imread('101087.jpg');
	spb = spbMIL(img);
	
For more details, take a look at the `spbDemo`.

*NOTE: Spectral feature extraction is time consuming. That is why we prefered to store the extracted features for the desired images and load them when necessary.*

## Citation 

If you find our code or annotations useful for your research, please cite our paper [Learning-Based Symmetry Detection for Natural Images](http://tsogkas.github.io/publications/symmetry-mil/tsogkas2012learning.pdf):

```
@inproceedings{tsogkas2012learning,
	title={Learning-based symmetry detection in natural images},
	author={Tsogkas, Stavros and Kokkinos, Iasonas},
	booktitle={European Conference on Computer Vision},
	pages={41--54},
	year={2012},
	organization={Springer}
}
```

## References

- Tsogkas, S., Kokkinos I.: Learning-Based Symmetry Detection in Natural Images. ECCV (2012)
- Martin, D., Fowlkes, C., Tal, D., Malik, J.: A database of human segmented natural
images and its application to evaluating segmentation algorithms and measuring
ecological statistics. In: ICCV (2001).
- Martin, D., Fowlkes, C., Malik, J.: Learning to detect natural image boundaries
using local brightness, color, and texture cues. PAMI (2004).
- Arbelaez, P., Maire, M., Fowlkes, C., Malik, J.: Contour detection and hierarchical
image segmentation. PAMI (2011)
- Telea, A., Van Wijk, J.: An augmented fast marching method for computing skeletons
and centerlines. In: Eurographics (2002).
- Kokkinos, I., Maragos, P., Yuille, A.: Bottom-up & top-down object detection using
primal sketch features and graphical models. In: CVPR, vol. 2, pp. 18931900.
IEEE (2006).

