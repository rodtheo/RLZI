# Relative Lempel-Ziv with Inexact Matchings (RLZI)

Briefly, it detects the substrings from a sample sequence S previously encountered in a reference sequence R with up to **k** mismatches (where k is a parameter). It belongs to the family of algorithms known as lempel-ziv parsers and its variants.

## Installing!

Clone this repository recursively to assure the download of submodules ([BWTIL](https://github.com/rodtheo/BWTIL) and [sdsl-lite](https://github.com/simongog/sdsl-lite)). Then, go into the repo directory and execute the install script. Furthermore, create the conda environment (`environment.yaml`) to install the necessary packages to execute the tests. Take a look at Quickstart section for detailed instructions.

### Quickstart (traditional):

Clone recursively the repository, enter the main folder, install RLZI, create the environment, activate it:

```
git clone --recursively https://github.com/rodtheo/RLZI
git submodule update --remote
cd ./RLZI/
bash install.sh
conda env create -f environment.yml -n rlzi-env
conda activate rlzi-env
```

### Quickstart (Docker):

```
git clone https://github.com/rodtheo/RLZI
cd ./RLZI/
docker run -it -v $PWD/test-data:/home/rlzi-user/RLZI/test-data rodtheo/rlzi:latest
(inside docker environment) conda activate /home/rlzi-user/env
(inside docker environment) snakemake --use-conda -p -j 1
```

### Optional (install from dockerfile)

```
git clone https://github.com/rodtheo/RLZI
cd ./RLZI/docker
docker build -t rlzi_img --file Dockerfile ../
```

## Running!

To run, execute (in the top-level directory):

```
snakemake --use-conda -p -j 1
```

This should succeed :).

MODIFY FOLLOWING SENTENCE.
_Once that works, you can configure it yourself by copying
`test-data/conf-test.yml` to a new file and editing it. See
`conf/conf-necator.yml` for a real example_

## Explanation of output files (still to be changed).

By default, RLZI will create the output directory in accordance to what is specified in the config file you use.

For instance, in `configs/config-test.yaml` are written the configurations to run the test dataset.

It will create the output directories `RLZ_k10_m1` and `RLZ_k10_m10` inside `test-data` because, in config file, we have declared the parameter `k=10` and parameter `m=10`. Therefore, RLZI will execute twice in order to create both output directories.

In the output directories (e.g. `test-data/RefL_TarL/RLZ_k10_m10`, or whatever is specified
in the config file you use), there will be a few important files --
the main ones are,

* `RLZ_k10_m10` -
* `RLZ_k10_m10_sdbB.sdsl` -
* `RLZ_k10_m10_sdbChref.sdsl` -
* `RLZ_k10_m10_sdbl.sdsl` -

## Resources

Parsing a 10 MB sample assembly in relation to a 10 MB reference assembly, RLZI took about 2 minutes and required about 116 MB of RAM.

## Need help?

Please ask questions and file issues on [the RLZI GitHub issue tracker](https://github.com/rodtheo/RLZI/issues).

## Credits

Thanks to ... for their inspiration!

This README is inspired in the files described at the excelent [blog](http://ivory.idyll.org/blog/2020-software-and-workflow-dev-practices.html) written by Dr. Titus Brown.

----
