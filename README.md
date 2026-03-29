# Shape Matching Benchmark Datasets

This repository provides a unified collection and preprocessing pipeline for widely used **shape matching benchmarks**, facilitating reproducible research and fair comparison across different methods.

---

## 📌 Overview

Shape matching aims to establish dense correspondences between 3D shapes under various transformations, including **non-rigid deformations**, **partiality**, **topological noise**, and **cross-category variations**. It is a fundamental problem in computer vision and computer graphics, with applications in:

- 3D reconstruction  
- Texture transfer  
- Animation and deformation modeling  
- Medical shape analysis  
- Cross-category semantic correspondence  

To support research in this area, we organize commonly used datasets into a **unified format**, enabling consistent evaluation and easier integration into existing pipelines.

---

## 📊 Supported Datasets

| Dataset | Category | Characteristics | Link |
|--------|----------|----------------|------|
| FAUST | Human | Near-isometric deformation | https://faust.is.tue.mpg.de |
| SCAPE | Human | Large non-rigid deformation | https://ai.stanford.edu/~drago/Projects/scape/scape.html |
| TOSCA | Synthetic | Canonical isometric benchmark | http://tosca.cs.technion.ac.il |
| Partial | Partial shape | Partial shape matching |https://cvg.cit.tum.de/data/datasets/partial |
| SMAL | Animal | Category-level variation | https://smal.is.tue.mpg.de |
| KIDS | Synthetic | Complex deformation | https://cvg.cit.tum.de/data/datasets/kids |
| TOPKIDS | Synthetic | Topological noise | https://cvg.cit.tum.de/data/datasets/topkids |
| SHREC16 | Partial matching | Partial shapes | https://www.shrec.net/ |
| SHREC19 | Cross-dataset | Generalization evaluation | https://www.shrec.net/ |
| SHREC20 | Partial / noisy | Realistic challenges | https://www.shrec.net/ |
| DT4D-H | Cross-category | Strong non-isometric deformation | https://github.com/rabbityl/DeformingThings4D |

---

## 📌 Citation

If you find this repository or the processed datasets helpful for your research, please consider citing our related works:

```bibtex
% Paper 1
@inproceedings{xia2024locality,
  title={Locality preserving refinement for shape matching with functional maps},
  author={Xia, Yifan and Lu, Yifan and Gao, Yuan and Ma, Jiayi},
  booktitle={Proceedings of the AAAI Conference on Artificial Intelligence},
  volume={38},
  number={6},
  pages={6207--6215},
  year={2024}
}

% Paper 2
@inproceedings{xia2025multi,
  title={Multi-shape matching with cycle consistency basis via functional maps},
  author={Xia, Yifan and Ye, Tianwei and Zhou, Huabing and Wang, Zhongyuan and Ma, Jiayi},
  booktitle={Proceedings of the AAAI Conference on Artificial Intelligence},
  volume={39},
  number={8},
  pages={8575--8583},
  year={2025}
}



% Paper 3
@article{xia2026locality,
  title={Locality Optimization Refinement with Deformation for Shape Matching via Functional Maps},
  author={Xia, Yifan and Ma, Jiayi},
  journal={International Journal of Computer Vision},
  volume={134},
  number={2},
  pages={76},
  year={2026},
  publisher={Springer}
}
