# SHADOW-TD v1.0
**Shadow Imaging with Thick Disks** 
is a Python package for general-relativistic ray-tracing imaging of black hole shadows with geometrically thick accretion disks.
It is lightweight and designed as an easy-to-use imaging tool that runs on standard personal computers, producing observer-plane images and numerical data for analyzing shadow structure and radiation flux distributions.
## Example Shadow Image
![Black Hole Shadow](https://private-user-images.githubusercontent.com/217470441/517449761-9d4cbb09-d211-4ba1-94c7-8e8f3935a3a1.png?jwt=eyJ0eXAiOiJKV1QiLCJhbGciOiJIUzI1NiJ9.eyJpc3MiOiJnaXRodWIuY29tIiwiYXVkIjoicmF3LmdpdGh1YnVzZXJjb250ZW50LmNvbSIsImtleSI6ImtleTUiLCJleHAiOjE3NjM3MzgzOTUsIm5iZiI6MTc2MzczODA5NSwicGF0aCI6Ii8yMTc0NzA0NDEvNTE3NDQ5NzYxLTlkNGNiYjA5LWQyMTEtNGJhMS05NGM3LThlOGYzOTM1YTNhMS5wbmc_WC1BbXotQWxnb3JpdGhtPUFXUzQtSE1BQy1TSEEyNTYmWC1BbXotQ3JlZGVudGlhbD1BS0lBVkNPRFlMU0E1M1BRSzRaQSUyRjIwMjUxMTIxJTJGdXMtZWFzdC0xJTJGczMlMkZhd3M0X3JlcXVlc3QmWC1BbXotRGF0ZT0yMDI1MTEyMVQxNTE0NTVaJlgtQW16LUV4cGlyZXM9MzAwJlgtQW16LVNpZ25hdHVyZT1iN2Q2Nzk3NDk4Zjc0NjMxZmI5NmEzMTZhNzE1YWJkNTQxNDJjYTUxNGRiN2I3NjU3ODZlOGFiNWQxODYwNzk5JlgtQW16LVNpZ25lZEhlYWRlcnM9aG9zdCJ9.9AxK3inKHuVsqjOVci8wEfkAdar_Qye5aTni50UxqLo)
## Features
1. Schwarzschild spacetime as the default background metric (to use a different spherically symmetric metric, modify the default metric in the `step1` and `step2` scripts within the `scripts/` folder).
2. Thick-disk illumination based on a generalized Shakura–Sunyaev model.
3. 2-D shadow images and flux distribution.
4. Editable configuration files in the `config/` folder for disk geometry, flow dynamics (rotation and infall velocities), opacity (absorption coefficient), and observer's viewing angle.

## Directory Structure
```
SHADOW-TD/
├── run_all.py # Main entry point
├── requirements.txt # Dependencies
├── README.md # Documentation
├── LICENSE # MIT License
├── config/ # Editable parameter files
├── scripts/ # Data processing and plotting
└── output/ # Images and numerical results
```
## Installation
For the first run, ensure that  all required packages are installed by executing the following command in the project's root directory:
 `pip install -r requirements.txt`

## Run the main pipeline 
`python run_all.py`
The program automatically loads configuration files from the `config/` directory, performs the ray-tracing computation, and saves all images and numerical data into `output/`.

## Custom Analysis
For flux distribution analysis or additional plots, run the corresponding scripts in `scripts/`.

## Citation
If you use this code in academic work, please cite:

Ziliang Wang, *Exploring the role of accretion disk geometry in shaping black hole shadows*,
Phys. Rev. D 112, 064052 (2025). arXiv:2506.21148[gr-qc]

## License
This project is released under the MIT License. See the `LICENSE` file for details.

## Author
Ziliang Wang,
Email: ziliang.wang@just.edu.cn
