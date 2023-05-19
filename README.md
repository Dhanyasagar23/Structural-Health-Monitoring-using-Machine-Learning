# Structural-Health-Monitoring-using-Machine-Learning
The focus of this project was on tackling inverse problems in the field of Structural Health Monitoring (SHM). The primary objective was to develop a methodology to estimate the stiffness of a simply supported truss structure comprising 20 members and 10 joints. 

To achieve this, a Finite Element Model (FEM) was constructed using MATLAB. The FEM allowed for the analysis of the structural behavior and the calculation of the stiffness based on the known nodal displacements. The FEM served as a reference model for evaluating the accuracy of the subsequent machine learning approach.

An Artificial Neural Network (ANN) was built as a machine learning model in MATLAB. The ANN was trained using the FEM data, and its purpose was to predict the stiffness of the truss members. By utilizing the ANN, it became possible to estimate the stiffness of the individual members without relying solely on the FEM calculations.

The trained ANN model was then employed to predict the stiffness of the truss members, and the results obtained were compared with the values obtained from the FEM. This comparison aimed to assess the performance and accuracy of the ANN model in predicting the stiffness values. 

The evaluation of the ANN model's performance was conducted using appropriate metrics and statistical measures to gauge the similarity between the predicted and FEM-derived stiffness values. The comparison helped determine the effectiveness of the ANN model as an alternative or supplementary method for estimating the stiffness in SHM applications.

