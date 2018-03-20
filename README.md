# MPI-MHD-FV

This is a numerical code developed to solve the magnetohydrodynamics equations. It was developed primarily for jets in the interstellar medium. 


### Status
[![Codacy Badge](https://api.codacy.com/project/badge/Grade/34b42d431caf402d8cb1c8bdb1dba221)](https://www.codacy.com/app/garethcmurphy/mpi-mhd-fv?utm_source=github.com&utm_medium=referral&utm_content=garethcmurphy/mpi-mhd-fv&utm_campaign=badger)
[![Build Status](https://travis-ci.org/garethcmurphy/mpi-mhd-fv.svg?branch=master)](https://travis-ci.org/garethcmurphy/mpi-mhd-fv)
[![DOI](https://zenodo.org/badge/50610732.svg)](https://zenodo.org/badge/latestdoi/50610732)





## Getting Started

These instructions will get you a copy of the project up and running on your local machine for development and testing purposes. See deployment for notes on how to deploy the project on a live system.

### Prerequisites

What things you need to install the software and how to install them

```
Docker
cmake (minimum 2.8.7)
gcc
ctest
Google Test
hdf5
```

### Installing

A step by step series of examples that tell you have to get a development env running

Say what the step will be

```
docker build -t mhd .
```

And repeat

```
docker run -it mhd /bin/bash 
```

You can use 
```
python3 ../see.py 
```
to plot a file

## Running the tests

Google Test/ctest platform is used

```
ctest -VV
```



## Deployment

Add additional notes about how to deploy this on a live system

## Built With


## Contributing

Please read [CONTRIBUTING.md](https://gist.github.com/PurpleBooth/b24679402957c63ec426) for details on our code of conduct, and the process for submitting pull requests to us.

## Versioning

We use [SemVer](http://semver.org/) for versioning. For the versions available, see the [tags on this repository](https://github.com/garethcmurphy/mpi-mhd-fv/tags). 

## Authors

* **Gareth Murphy** - *Initial work* - [garethcmurphy](https://github.com/garethcmurphy)

See also the list of [contributors](https://github.com/your/project/contributors) who participated in this project.

## License

This project is licensed under the MIT License - see the [LICENSE.md](LICENSE.md) file for details

## Acknowledgments

* Hat tip to anyone who's code was used
* Inspiration
* etc
