# New methods for bird sister species comparison

QDR, 3 Feb 2017

## Abstract

In 1967, Dan Janzen wrote a clever paper that was far ahead of its time. He stated that mountain passes are (physiologically speaking) higher in the tropics because the lower temporal variation experienced by tropical species causes them to be narrowly adapted to local conditions and thus have narrower niches. This idea has some support when it comes to thermal niches, but it has rarely been tested for niches in general. Body size is a good proxy for many niche dimensions in vertebrates. We identified sister species pairs using a phylogeny of all the birds in the world. For each of these pairs, we searched the VertNet specimen database for pairs where each species had at least ten body mass measurements. We tested the hypothesis that tropical species have less intraspecific variability in body mass--indeed, it is supported. This finding provides strong evidence for an important macroecological hypothesis. In addition, the lower standing variation in body mass of tropical species hints at narrower niche breadths, which may make them less able to cope with changing environments in the future. 

## New methods

Instead of beginning with the VertNet database, this time I began with the phylogeny.

### Main comparison

1. I downloaded the phylogeny of all birds from birdtree.org and took a random sample of 50 phylogenetic trees.
2. I calculated the pairwise branch length between all pairs of taxa in each of the trees, then averaged the 50 distance matrices.
3. I defined a sister species pair as a pair of species A and B, such that both the following are true: species A is species B's single nearest neighbor, and species B is species A's single nearest neighbor. This should represent two species that are each other's closest phylogenetic relative. I found all such pairs in the phylogeny (approximately 2700 pairs).
4. I downloaded all bird body mass records from vertnet.org. I ran some basic quality control on the records, removing anomalous records and all non-adult birds.
5. I excluded all pairs in which any species had less than ten valid body mass measurements (leaving approximately 700 pairs).
6. I used the median latitude of all the records for each species to classify each species as either tropical (`$|lat| < 23.5$`) or temperate (`$|lat| > 23.5$`). I excluded all pairs except for tropical-temperate pairs. This left 101 pairs.
7. I calculated the coefficient of variation (CV) of `$\log_{10} mass$` for each species.
8. I ran a paired, one-sided t-test on the species pairs. The alternative hypothesis, based on theory, is that the CV of the tropical species is lower.
9. To make sure that there was no effect of differing sample sizes within species pairs, I used a bootstrap procedure as follows: I took a subsample, without replacement, of the species with more records in each pair. The subsample size was equal to the number of records of the other species in the pair. I recalculated the CV of each species and ran the t-test. I repeated this procedure 999 times to get a bootstrap distribution of t-statistics and p-values. 

### Covariates

We were additionally interested in testing the influence of different covariates on the pattern.

1. *Higher taxonomic groupings*: 62 of the 101 pairs were passerines, and the remaining 39 were in a grab bag of orders (too few per family/order to say anything meaningful about individual orders). I grouped the pairs by passerine/non-passerine and tested whether the mean difference was different between the two groups.
2. *Functional traits*: I extracted foraging traits from a database of bird functional traits (EltonTraits, http://dx.doi.org/10.1890/13-1917.1). I extracted life-history traits from a database of amniote life history traits (https://doi.org/10.6084/m9.figshare.c.3308127.v1). I tested whether the mean difference varied by the literature value of species mean body mass, or by a categorical diet classification. 
3. *Range where the specimens were collected*: I used the lat-long coordinates of each point where specimens were collected to calculate a very crude measurement of the spatial spread of the collected specimens (the mean pairwise distance between the coordinates in lat-long space, this can be refined later). I tested whether the mean difference was correlated with this spread.
4. *Climate variability where the specimens were collected*: I downloaded the global 1979-2013 BioClim variables for the entire terrestrial globe at 0.5' resolution (http://dx.doi.org/doi:10.1594/WDCC/CHELSA_v1_1). I calculated the following variables and tested whether the mean difference was correlated with them: 
	- spatial CV of temperature across the range where specimens were collected
	- spatial CV of precipitation across the range
	- mean seasonal variability of temperature, averaging all the points where specimens were collected (temporal var within year)
	- mean seasonal variability of precipitation (temporal var within year)
	- mean interannual variability of temperature (temporal var across years)
	- mean interannual variability of precipitation (temporal var across years)
5. *Range sizes*: I have downloaded range sizes from BirdLife International at birdlife.org, but I am still working on extracting useful information out of them.
6. *Migrant status*: I have some information from the Clements Checklist, but it is incomplete in terms of species. I might be able to get it off of the range maps. This is still a work in progress as well.

## Results

Tropical bird species have lower CV of body mass than their temperate sister species: `$t_{100} = -3.06$`, `$p = 0.0015$`, mean difference `$\bar {d} = -0.020$` (less in the tropics), `$CI = (-\infty, -0.0092)$`. (Figure 1: Plot connecting each pair with a line; Figure 2: Histogram of differences).

This result held up to the bootstrapping test, indicating that it is not an artifact of different sample sizes between tropical and temperate birds. The t-test indicated that the tropical species had a significantly lower CV in 95.1% of the bootstrap simulations. The central 95% of the distribution of mean differences from the simulations was `$(-0.024, -0.010)$` and the central 95% of the distribution of t-statistics was `$(-3.22, -1.53)$`.

## To do

I will compile information on the functional group (using published bird trait databases), range size (using birdlife's range maps), and migrant status (using the eBird Clements checklist, as well as other sources as needed). Each of these can be used to look at whether the pattern between tropical and temperate is mediated by functional group, range size, or migrant status. I will also get the location of each specimen and look at whether the CV of a species is a function of the variability in collection location. I also can get climate data for all the specimen locations. How should that be put into this?

<div style="page-break-after: always;"></div>
## Figures

![bird figure](file:///C:\\Users\\Q\\Dropbox\\projects\\verts\\vertnet_results\\birdinteractionplot03feb.png)  
**Figure 1.** Interaction plot showing that the tropical species has lower intraspecific variability in body mass than the non-tropical species. Each point is the CV of log-transformed body mass for each species. Sister taxa are connected with a gray line. The red points are the mean CVs for each realm +/- 1 SE, connected by a red line. 

![histogram](file:///C:\\Users\\Q\\Dropbox\\projects\\verts\\vertnet_results\\birdhistogram03feb.png) 
**Figure 2.** Histogram of the mean differences for birds showing they are normally distributed and that the paired t-test is appropriate. The thick red line is the mean difference, the thin red line is the upper bound of the confidence interval, and the dotted blue line is zero. This shows that the mean difference is significantly less than zero.