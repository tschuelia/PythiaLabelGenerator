# Pythia Difficulty Label Generator

> [!CAUTION]
> 🚧 **Under Construction** 🚧
> This project is still in development and is not yet ready for use.

The Pythia Difficulty Label Generator is a tool that generates the ground-truth phylogenetic difficulty label for an
MSA.

> [!WARNING]
> Computing the ground-truth difficulty for an MSA is very time-consuming and requires a lot of computational resources,
> especially for MSAs with many sites and/or taxa.
> Please use our (very accurate) difficulty prediction tool [Pythia](https://github.com/tschuelia/PyPythia) instead
> whenever possible.

## Difficulty Label Computation

This tool computes the difficulty according to our definition published in [Haag _et
al._ (2022)](https://doi.org/10.1093/molbev/msac254).

Let $N_{\text{all}}$ be the number of inferred ML trees.
We first compute the average pairwise relative Robinson-Foulds (RF) Distance between all trees ($RF_{\text{all}}$), as
well as the number of unique tree topologies among the inferred trees ($N^*_{\text{all}}$).
We filter the inferred trees using likelihood-based statistical tests to obtain the set of $N_{\text{pl}}$ _plausible
trees_.
We again compute the average pairwise RF distance between the plausible trees ($RF_{\text{pl}}$) and the number of
unique tree topologies among the plausible trees ($N^*_{\text{pl}}$).

The difficulty is then computed as follows:

$$
\text{difficulty} = \frac{1}{5} \cdot \bigg[ RF_{\text{all}} + RF_{\text{pl}}
+ \frac{N^*_{\text{all}}}{N_{\text{all}}} + \frac{N^*_{\text{pl}}}{N_{\text{pl}}}
+ \left( 1 - \frac{N_{\text{pl}}}{N_{\text{all}}} \right) \bigg]
$$

For further details on the reasoning and validation of this difficulty, please refer to the publication linked above.

We infer the ML trees using [RAxML-NG]((https://github.com/amkozlov/raxml-ng) and use the statistical significance tests
as implemented in [IQ-TREE](http://www.iqtree.org).
Per default, the difficulty is based on $N_{\text{all}}=100$ ML trees.
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

