# Link Prediction and Navigability of Multiplex Energy Networks

## Overview
This project applies advanced link prediction techniques to multiplex networks, specifically focusing on energy networks. It features an innovative approach that enhances traditional link prediction algorithms for better accuracy in multi-layered networks. The codebase includes simulations based on Belgium's electricity and gas networks, demonstrating how different random walk strategies and strategic link additions can influence network navigability and resilience. This project is ideal for researchers and practitioners interested in network analysis, link prediction, and network resilience strategies.

### Codes
- **Link_Prediction.ipynb**: This Jupyter notebook contains the code for predicting new edges in multiplex networks. It includes the algorithms and methods necessary for link prediction, providing a comprehensive guide to the process.
- **Navigability_RW_Strategies.R**: This script contains various navigability and random walk strategies that can be applied to the predicted edges. While the project currently implements classical, PageRank, and diffusive strategies, the code is also equipped to handle other random walk strategies, depending on the feasibility of the network.

### Data
The data directory contains five Excel sheets, each representing a different layer of the multiplex network. These sheets are crucial for the analysis and are used by the code to perform link prediction and navigability studies.

- **Case_1.xlsx**
- **Case_2.xlsx**
- **Case_3.xlsx**
- **Case_4.xlsx**
- **Case_5.xlsx**

## Getting Started
To use this repository:
1. Clone the repository to your local machine.
2. Ensure you have the necessary dependencies installed, including Python for the Jupyter notebook and R for the R scripts.
3. Start with the `Link_Prediction.ipynb` notebook to understand the link prediction process.
4. Explore the `Navigability_RW_Strategies.R` script to see how different navigability strategies can be applied to your network.
5. Use the data provided in the Data directory to test and apply the code.

## Contributing
Contributions to this project are welcome. Please feel free to fork the repository, make your changes, and submit a pull request.

## Citation
If you use the code or data in this repository in your work, please cite the following paper:

Kazim, M., Pirim, H. (2024). Multilayer Analysis of Energy Networks. *ASC*. Available at: 

## Contact
muhammad.kazim@ndsu.edu
