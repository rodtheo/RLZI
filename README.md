# Relative Lempel-Ziv with Inexact Matchings (RLZI)

Briefly, it detects the substrings from a sample sequence S previously encountered in a reference sequence R with up to $k$ mismatches (parameter). It belongs to the family of algorithms known as lempel-ziv parsers and its variants.

## Installing!

Clone this repository recursively to assure the download of submodules ([BWTIL](https://github.com/rodtheo/BWTIL) and [sdsl-lite](https://github.com/simongog/sdsl-lite)). Then, go into the repo directory and execute the install script. Furthermore, create the conda environment (`environment.yaml`) to install the necessary packages to execute the tests. Take a look at Quickstart section for detailed instructions.

### Quickstart:

Clone recursively the repository, enter the main folder, install RLZI, create the environment, activate it:

```
git clone --recursively https://github.com/rodtheo/RLZI
cd ./RLZI/
bash install.sh
conda env create -f environment.yml -n rlzi-env
conda activate rlzi-env
```

## Running!

To run, execute (in the top-level directory):

```
snakemake --use-conda -p -j 1
```

This should succeed :).

Once that works, you can configure it yourself by copying
`test-data/conf-test.yml` to a new file and editing it. See
`conf/conf-necator.yml` for a real example.

## Explanation of output files (still to be changed).

In the output directory (e.g. `output.test`, or whatever is specified
in the config file you use), there will be a few important files --
the main ones are,

* `gather.csv` - the list of contaminants
* `matching-contigs.fa` - all contigs with any matches to the database
* `matching-fragments.fa` - all fragments with any matches to the database

## Resources

Parsing a 10 MB sample assembly in relation to a 10 MB reference assembly, RLZI took about 2 minutes and required about 116 MB of RAM.

## Need help?

Please ask questions and file issues on [the RLZI GitHub issue tracker](https://github.com/rodtheo/RLZI/issues).

## Credits

Thanks to ... for their inspiration!

This README is inspired in the files described at the excelent [blog](http://ivory.idyll.org/blog/2020-improved-workflows-as-applications.html) written by Dr. Titus Brown.

----
