# Relative Lempel-Ziv with Inexact Matchings (RLZI)

Briefly, it detects the substrings from a sample sequence S previously encountered in a reference sequence R with up to **m** mismatches spaced by at least **k** characters (where k and m are parameters). It belongs to the family of algorithms known as lempel-ziv parsers and its variants.

## Installing!

### Recommended (using Docker):

Make sure you have [Docker](https://www.docker.com/) and [Conda](https://docs.anaconda.com/anaconda/install/) (or [Miniconda3](https://docs.conda.io/en/latest/miniconda.html)) installed and running smoothly in your computer.

After this, clone the repository, enter the main folder, create and activate the python3.x environment:

```
git clone https://github.com/rodtheo/RLZI
cd ./RLZI/
conda env create -f environment.yml -n rlzi-env
docker pull rodtheo/rlzi:v1
```

### Install from source (compiling the software):

#### Dependencies
Boost
To add more..

Clone this repository recursively to assure the download of submodules ([BWTIL](https://github.com/rodtheo/BWTIL) and [sdsl-lite](https://github.com/simongog/sdsl-lite)). Then, go into the repo directory and execute the install script. Furthermore, create the conda environment (`environment.yaml`) to install the necessary packages to execute the tests. 

```
git clone --recursively https://github.com/rodtheo/RLZI
git submodule update --remote
cd ./RLZI/
bash install.sh
conda env create -f environment.yml -n rlzi-env
```

## Running!

To run, activate the environment `rlzi-env` and execute the `Snakemake` file (in the top-level directory):

```
conda activate rlzi-env
snakemake --use-conda -p -j 1
```

This should succeed :).

Once that works, you can configure it yourself by copying `configs/config-test.yaml` to a new file and editing it.

## Explanation of output files (_still to be changed_).

By default, RLZI will create the output directory in accordance to what is specified in the config file you use.

For instance, in `configs/config-test.yaml` are written the configurations to run the test dataset.

It will create the output directories `RLZ_k10_m1` and `RLZ_k10_m10` inside `test-data/RefS_TarL` because, in config file, we have declared the parameter `k=10` and parameters `m={1, 10}`. Therefore, RLZI will execute twice in order to create both output directories.

The created folders will be `test-data/RefL_TarL/RLZ_k10_m1` and `test-data/RefL_TarL/RLZ_k10_m10`. Let's examine the files within one of the output dirs.

For instance, in `test-data/RefL_TarL/RLZ_k10_m10` (or whatever is specified in the config file you use), there will be a few important files, the main ones are:

* `RLZ_k10_m10` -
* `RLZ_k10_m10` - 

## Resources

Parsing a 10 MB sample assembly in relation to a 10 MB reference assembly, RLZI took about 2 minutes and required about 116 MB of RAM.

## Need help?

Please ask questions and file issues on [the RLZI GitHub issue tracker](https://github.com/rodtheo/RLZI/issues).

## Credits

Thanks to ... for their inspiration!

This README is inspired in the files described at the excelent [blog](http://ivory.idyll.org/blog/2020-software-and-workflow-dev-practices.html) written by Dr. Titus Brown.

----
