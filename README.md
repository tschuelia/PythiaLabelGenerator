# Pythia Difficulty Label Generator

![Label Generator GH actions CI](https://github.com/tschuelia/PythiaLabelGenerator/actions/workflows/test-label-generator.yml/badge.svg)

The Pythia Difficulty Label Generator generates the ground-truth phylogenetic difficulty label for an MSA and
corresponds to the prediction target of our difficulty prediction tool [Pythia](https://github.com/tschuelia/PyPythia).

> [!CAUTION]
> Computing the ground-truth difficulty for an MSA is very time-consuming and requires a lot of computational resources,
> especially for MSAs with many sites and/or taxa.
>
> Please use our (very accurate) difficulty prediction tool [Pythia](https://github.com/tschuelia/PyPythia) instead
> whenever possible.

The ground-truth difficulty is computed according to our definition published in [Haag _et
al._ (2022)](https://doi.org/10.1093/molbev/msac254):

Let $`N_{\text{all}}`$ be the number of inferred ML trees.
We first compute the average pairwise relative Robinson-Foulds (RF) distance between all trees ($`RF_{\text{all}}`$), as
well as the number of unique tree topologies among the inferred trees ($`N^*_{\text{all}}`$).
We filter the inferred trees using likelihood-based statistical tests to obtain the set of $`N_{\text{pl}}`$ _plausible
trees_.
We again compute the average pairwise RF distance between the plausible trees ($`RF_{\text{pl}}`$) and the number of
unique tree topologies among the plausible trees ($`N^*_{\text{pl}}`$).

The difficulty is then computed as follows:

```math
\text{difficulty} = \frac{1}{5} \cdot \bigg[ RF_{\text{all}} + RF_{\text{pl}}
+ \frac{N^*_{\text{all}}}{N_{\text{all}}} + \frac{N^*_{\text{pl}}}{N_{\text{pl}}}
+ \left( 1 - \frac{N_{\text{pl}}}{N_{\text{all}}} \right) \bigg]
```

For further details on the reasoning and validation of this difficulty, please refer to the publication linked above.

We infer the ML trees using [RAxML-NG](https://github.com/amkozlov/raxml-ng) and use the statistical significance tests
as implemented in [IQ-TREE](http://www.iqtree.org).
Per default, the difficulty is based on $`N_{\text{all}}=100`$ ML trees.
Note that this number can be adjusted by the user, however, the difficulty will only be an approximation if the number
of trees is changed.

## Prediction of Phylogenetic Difficulty

As stated above, computing the ground-truth difficulty for an MSA is very time-consuming and requires a lot of
computational resources, especially for MSAs with many sites and/or taxa.
Please use our (very accurate) difficulty prediction tool [Pythia](https://github.com/tschuelia/PyPythia) instead
whenever possible.
For instance, inferring _a single ML tree_ for an MSA
comprising [SARS-CoV-2 sequences](https://doi.org/10.1093/molbev/msaa314) (approx. 5k taxa and 28.5k sites) takes about
12 hours on a large compute cluster. Using Pythia instead, we can predict the same MSA to be very difficult in about 2.5
minutes on a standard MacBook.

Only use this tool if you need the ground-truth difficulty for a specific MSA and you are sure that Pythia is unable to
predict the difficulty accurately.
The only case where we observed Pythia to fail is for language MSAs, so if you are working with DNA, Protein, or
biological morphological data, Pythia should work just fine 😉

## Installation

#### Requirements

To use this labelling tool, you need to install

- RAxML-NG: See [the RAxML-NG GitHub repository](https://github.com/amkozlov/raxml-ng) for installation instructions.
  Please make sure that you install a RAxML-NG version < 2.
- IQ-TREE: See [the IQ-TREE website](http://www.iqtree.org) for installation instructions. Please install IQ-TREE
  version 2 or higher.

#### Install via conda (recommended)

This package will soon be available on conda-forge :)

#### Install using pip

You can install the package using pip:

```bash
pip install pythialabelgenerator
```

## Usage

This label-generator is primarily a command line tool. You can call it using the `label` command, for instance, to
compute the difficulty for the example MSA provided in the `examples` directory run

```bash
label -m examples/example.phy
```

This will infer 100 ML trees using RAxML-NG and run the statistical tests using IQ-TREE. The difficulty will be printed
to the console.
The output will look something like this:

```text
Difficulty LabelGenerator version 1.0.0 released by The Exelixis Lab
Developed by: Julia Haag
Latest version: https://github.com/tschuelia/LabelGenerator
Questions/problems/suggestions? Please open an issue on GitHub.

LabelGenerator was called at 06-Mar-2025 15:15:03 as follows:

label -m examples/example.phy

[00:00:00] Starting label computation.
[00:00:00] Inferring 100 ML trees using RAxML-NG.
[00:00:23] Computing RF-Distance between ML trees.
[00:00:23] > RF-Distance ML trees: 0.78
[00:00:23] > Unique topologies ML trees: 100
[00:00:23] Running IQ-TREE statistical tests.
[00:00:27] Filtering plausible ML trees.
[00:00:27] > Found 18 plausible trees.
[00:00:27] Computing RF-Distance between plausible ML trees.
[00:00:27] > RF-Distance plausible trees: 0.78
[00:00:27] > Unique topologies plausible trees: 18

Ground Truth Difficulty for examples/example.phy: 0.875

Total runtime: 27.42 seconds.
```

Depending on your system setup, you might need to pass a RAxML-NG and IQ-TREE binary path to the label generator.
You can do this using the `-r` and `-i` options, respectively. This is required in case `raxml-ng` and/or `iqtree2` are
not in your `$PATH`.

Note that this `examply.phy` MSA is not the same exemplary MSA as we provide in the PyPythia repository, so please don't compare this ground-truth lable to the exemplary prediction in PyPythia 😉

For a full list of command line options, run `label -h`:

```text
Difficulty LabelGenerator version 1.0.0 released by The Exelixis Lab
Developed by: Julia Haag
Latest version: https://github.com/tschuelia/LabelGenerator
Questions/problems/suggestions? Please open an issue on GitHub.

usage: label [-h] -m MSA -r RAXMLNG -i IQTREE [-t THREADS] [-s SEED] [-p PREFIX] [--model MODEL] [--ntrees NTREES] [--redo] [-V]

Generate the ground truth difficulty for the given MSA.

options:
  -h, --help            show this help message and exit
  -m MSA, --msa MSA     Multiple Sequence Alignment to compute the ground truth difficulty for. Must be in either phylip or fasta format.
  -r RAXMLNG, --raxmlng RAXMLNG
                        Path to the binary of RAxML-NG. For install instructions see https://github.com/amkozlov/raxml-ng.(default: 'raxml-
                        ng' if in $PATH, otherwise this option is mandatory).
  -i IQTREE, --iqtree IQTREE
                        Path to the binary of IQ-TREE2. For install instructions see http://www.iqtree.org.(default: 'iqtree2' if in $PATH,
                        otherwise this option is mandatory).
  -t THREADS, --threads THREADS
                        Number of threads to use for the RAxML-NG tree inference and IQ-TREE statistical tests (default: autoconfig in RAxML-NG and IQ-TREE).
  -s SEED, --seed SEED  Seed for the RAxML-NG tree inference (default: 0).
  -p PREFIX, --prefix PREFIX
                        Prefix of the RAxML-NG and IQ-TREE log and result files (default: MSA file name).
  --model MODEL         Model to use for the RAxML-NG tree inference (default: 'GTR+G' for DNA, 'LG+G' for AA, and 'MULTIx_GTR' for
                        morphological data where x is the maximum state value in the MSA).
  --ntrees NTREES       Number of ML trees to infer (default: 100)
  --redo                Redo all computations, even if the results already exist.
  -V, --version         Print the version number and exit.
```

Please note that inferring 100 ML trees for a large MSA can take a long time. You can adjust the number of trees to
infer using the `--ntrees` option. However, using fewer than 100 trees will likely result in slight different
difficulties, as the difficulty is based on the average pairwise RF distance between the inferred trees.

### Result Files

Running this labelling tool will result in the following files:

- `{prefix}.raxml.*`: RAxML-NG log and result files.
- `{prefix}.iqtree.*`: IQ-TREE log and result files.
- `{prefix}.labelGen.log`: Log file containing the output of the label generator. This is the same output as printed to
  the terminal.

You can set the prefix of these files using the `-p` option. By default, the prefix is the name of the MSA file.
Note that RAxML-NG and IQ-TREE refuse to overwrite existing files. If you want to redo the computations, you can use the
`--redo` option.
Please also specify the `--redo` option if you want to change the number of trees to infer using the `--ntrees` option
for the same prefix. Otherwise,
the label generator will exit with an error message.

### Input Data

You can provide the MSA in either phylip or fasta format. We currently support DNA, AA, and categorical data.
The label generator will automatically determine the data type and select an appropriate model for tree inference and
statistical tests.
For DNA data, we use the `GTR+G` model, for AA data the `LG+G` model, and for categorical data the `MULTIx_GTR` model, where `x`
is the maximum state value in the MSA.
We used these models to generate the ground-truth labels for training Pythia. If you want to specify a different model,
you can do so using the `--model` option.


## Citation
We will soon publish a pre-print on bioRxiv with updates on our Pythia difficulty prediction tool that will also include a
brief description of this new labelling tool. Please cite this pre-print if you use this tool in your research.

The link to the paper will be added soon 🙂
