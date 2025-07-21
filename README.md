# Code_Case5
This is the main computational code of the paper
"Determining the completed bridge state of a self-anchored suspension bridge: A non-iterative analytical approach"

# Project Description
This paper proposes a novel method based on design parameters such as 
control point coordinates (anchor points and IP points), sag-to-span ratios, and hanger spacing
to establish expressions for internal forces and deformations throughout the entire bridge using
a set of fundamental unknowns. By simultaneously considering closure conditions, static equilibrium conditions,
and zero vertical deflection constraints at hanger points, a closed system of nonlinear equations strictly 
matching the number of fundamental unknowns is formulated. Using the code in this project to solve this nonlinear equation system can obtain the fundamental unknowns.

# Instructions
test1 contains all files for Case 5 in the paper, with the following file purposes:

1. Ansys_Modeling_coordinates_01m - Modeling data for ANSYS 2022 R1:
(1) Cor_Grid_beam_01m.txt: Modeling coordinates for the main girder
(2) Cor_hanger_B.txt: Modeling coordinates for downstream hangers
(3) Cor_hanger_F.txt: Modeling coordinates for upstream hangers
(4) Cor_ZhuLan_B.txt: Modeling coordinates for upstream main cable
(5) Cor_ZhuLan_F.txt: Modeling coordinates for downstream main cable
(6) Real_Beam_01m.txt: Real constants for main girder
(7) Real_Steel.txt: Real constants for hangers
(8) Real_ZhuLan.txt: Real constants for main cable
(9) yugu_shang.txt: Modeling coordinates for upstream fishbone truss
(10) yugu_xia.txt: Modeling coordinates for downstream fishbone truss
2. Case5.m: Main computational code
3. chuzhi_3D_hangerXLX.txt: Initial values for the main code
4. Input files required for main code computation:
5. Coef_Eta4.txt: Coefficient data of the stiffening girder's design geometric shape
6. qg.txt: Distributed load parameters
7. xP.txt: Coordinate parameters

# Contact Information
You can achieve this by contacting 230239453@seu.edu.cn.