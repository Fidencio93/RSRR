# Design and Computational Modeling of a Leg-Wheel Transformable Mechanism with Decoupled Kinematics

**Authors**
* José F. Flores
* Héctor A. Moreno
* Isela G. Carrera
* José Luis Ordoñez

## Scripts
This Repository contains all the scripts required for the simulation and numerical results presented in the article

| Script | Description | Required Files | Related figures or tables |
| --- | --- | --- | --- |
| `NR.m` | Computes error using Newton-Raphson algorithm to try and solve inverse kinematics | `Test.mat` | Figure 4 |
| `Error_NN.m` | Calculates error using a Neural network to solve inverse kinematics   | `Cin_dir_2.m`, `Fcn_Red_Qp_1.m`, `Fcn_Red_Qr_1.m`,`Fcn_Red_Qp_2.m`, `Fcn_Red_Qr_2.m`,`Fcn_Red_Qp_3.m`, `Fcn_Red_Qr_3.m`,`Fcn_Red_Qp_4.m`, `Fcn_Red_Qr_4.m`,`Fcn_Red_Qp_5.m`, `Fcn_Red_Qr_5.m`,`Fcn_Red_Qp_6.m`, `Fcn_Red_Qr_6.m`, `Test.mat`  | Figure 5
| `Simulacion_cinematica.m` | Simulates the movement of the mechanism following a trajectory and plotting _x_ and _y_ coordinates| `Fcn_Red_Qp_1.m`, `Fcn_Red_Qr_1.m`, `Walk_50.mat` | Figure 7 |
| `Simulacion_cinematica_1.m` | Simulates the movement of the mechanism following a trajectory and plotting _x_ and _y_ coordinates| `Fcn_Red_Qp_1.m`, `Fcn_Red_Qr_1.m`, `Walk_50.mat` | Figure 6 |
| `Simulacion_cinematica_2.m` | Simulates the movement of the mechanism following a trajectory and plotting _x_ and _y_ coordinates| `Fcn_Red_Qp_1.m`, `Fcn_Red_Qr_1.m`, `Walk_50.mat` | Figure 6 |
| `M_perf.py` | Measures average execution time of neural networks |`Funcion_Qp_1.py`,`Funcion_Qr_1.py`, `Fcn_Qp.mat`,`Fcn_Qr.mat`,`Walk_50.mat` | NA |

## Requirements for matlab scripts
* Matlab 2018r or  later.
* No additional toolboxes are required.

## Requirements for python scripts
* Raspberry pi 4B (Only to measure NN performance)
* Python 3
*  numpy
*  scipy.io
*  time
*  statistics
