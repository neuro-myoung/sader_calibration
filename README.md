# Sader Calibration of Cantilever Stiffness

## Description
This is a python package written to calibrate cantilevers, for use in atomic force microscopy, using the method of Sader (Sader et al., 1998; Sader et al., 1999; Chon et al., 2000). This program uses the nidaqmx package to interface with an NI USB-6361 Data Aquisition Board (National Instruments). 

## Prerequisites

Before you begin make sure you have Python 3.9 or higher. It may work with other versions but I have not tested it with these particular packages.

## Installation

```
git clone https://github.com/neuro-myoung/sader_calibration.git
```

The environment.yml file has all the necessary packages and appropriate versions. Using conda create a virtual environment with the necessary packages using the environment.yml file by going to the folder in your command line and typing the following command:

```
conda env create -f environment. yml
```

The environment will be named saderCalibration by default. Activate the new environment whenever you want to run the program by entering in your command line the following:

```
conda activate saderCalibration
```

With the environment active you should be able to run the program by typing the following command from inside the package folder:

```
python sader_calibration.py
```

## Demo

First, input the number of samples per spectra, the total number of spectra to acquire, and the sampling rate. This program will acquire multiple spectra (I use 512 as a starting point) and use the method of Welch to average the spectra using hamming windows in order to reduce noise (Welch 1967). Once these are in click the **Acquire** button and you should see the analog inputs in the top plot and the continuously updating mean power spectra in the lower plot as seen in the gif below.

![Demo GIF](/assets/sampling.gif)

Once acquisition has finished, input the plane view dimensions of the cantilever (length and width) into the appropriate input boxes. Click the **Select** button which will now allow you to drag over the region of the plot corresponding to the resonance peak of the cantilever. Next, click the **Fit** button and it will fit a single harmonic oscillator (SHO) equation to the resonance peak. The table on the bottom left should then update with the fit parameters. Finally, click the **Calculate** button and it will determine the cantilever stiffness from the fit parameters and the plane view dimensions using the method of Sader (Sader et al., 1998; Sader et al., 1999; Chon et al., 2000). This second part is demonstrated in the gif below.

![Demo GIF2](/assets/calculation.gif)

## Contributing
To contribute to **sader_calibration**, follow these steps:

1. Fork this repository.
2. Create a branch: git checkout -b *branch_name*.
3. Make your changes and commit them: git commit -m '*commit_message*'
4. Push to the original branch: git push origin *project_name* *location*
5. Create a pull request.

Alternatively see the GitHub [documentation](https://help.github.com/en/github/collaborating-with-issues-and-pull-requests/creating-a-pull-request) on creating a pull request.

## Contributors

[@neuro-myoung](https://github.com/neuro-myoung)

## Contact

If you want to contact me you can reach me at michael.young@duke.edu

## License
This project uses an [MIT License](https://opensource.org/licenses/MIT)