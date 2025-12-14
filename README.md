# SHADOW-TD v1.0
**Shadow Imaging with Thick Disks** 
is a Python package for general-relativistic ray-tracing imaging of black hole shadows with geometrically thick accretion disks.
It is lightweight and designed as an easy-to-use imaging tool that runs on standard personal computers, producing observer-plane images and numerical data for analyzing shadow structure and radiation flux distributions.
## Example Shadow Image
![Black Hole Shadow](https://github.com/user-attachments/assets/78766ae0-d20a-4186-8dbb-02f9090a10e0)

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



