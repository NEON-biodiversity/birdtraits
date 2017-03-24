# New methods for bird sister species comparison

QDR   
created 3 Feb 2017  
modified 7 Feb 2017 (added info on covariates)  
modified 23 Mar 2017 (expanded discussion of covariates)

## Summary

In 1967, Dan Janzen wrote a clever paper that was far ahead of its time. He stated that mountain passes are (physiologically speaking) higher in the tropics because the lower temporal variation experienced by tropical species causes them to be narrowly adapted to local conditions and thus have narrower niches. This idea has some support when it comes to thermal niches, but it has rarely been tested for niches in general. Body size is a good proxy for many niche dimensions in vertebrates. We identified sister species pairs using a phylogeny of all the birds in the world. For each of these pairs, we searched the VertNet specimen database for pairs where each species had at least ten body mass measurements. We tested the hypothesis that tropical species have less intraspecific variability in body mass&mdash;indeed, it is supported. This finding provides strong evidence for an important macroecological hypothesis. In addition, the lower standing variation in body mass of tropical species hints at narrower niche breadths, which may make them less able to cope with changing environments in the future. 

The pattern of reduced variability in tropical taxa demands a mechanistic explanation. Janzen's original explanation rests on reduced intra-annual (seasonal) temperature variation in the tropics relative to temperate areas. Since an organism in the lowland tropics experiences a small range of temperatures within one year, there is no selective pressure to adapt to a wide range of potential temperatures. We argue that there is no reason to limit this hypothesis to within-year variation: tropical climates also tend to have lower interannual variation, which could also contribute to reduced niche breadth in tropical species. Our results provide moderate support (significant but relatively low explanatory power) for the hypothesis that reduced niche breadth in tropical species is driven by reduced climate variability: temperate-tropical bird pairs that have a higher discrepancy in temperature seasonality and in interannual temperature variation also tend to have a higher discrepancy in body mass variability.

Other mechanisms could explain the difference between temperate and tropical bird body mass variability. For instance, Rapoport's Rule states that tropical species tend to have smaller ranges; this could be correlated with reduced size variation. However, we found no influence of range size. In addition, classical niche theory predicts that coexistence of many species is permitted by individual species reducing their niche breadths to "pack" more species into niche space. To test this, we also looked at the effect of overall bird species richness, and of the number of coexisting congeners (both roughly calculated from global species range maps), but found no relationship in either case. Intriguingly, we found that there were greater relative temperate-tropical differences in large-bodied species; again, this is a significant pattern with relatively low explanatory power. 

We also included a number of other covariates that are essentially artifacts of data collection. These included the size of the geographic area from which the museum specimens were collected, the spatial variability in climate over which the specimens were collected, the number of populations from which the specimens were collected, and the elevational range from which they were collected. None of these potentially confounding factors had any relationship with the difference in body mass variability among temperate-tropical species pairs. (This is good because those patterns wouldn't be interesting!) 

**Take-home message: Tropical birds have a narrower range of body sizes than temperate sister taxa. This supports Janzen's hypothesis. Furthermore, the only plausible mechanism that explains any of the variation is differences in temperature variability (we cannot tell whether interannual or seasonal is more important).**

<div style="page-break-after: always;"></div>
## Methods

Instead of beginning with the VertNet database, this time I began with the phylogeny.

### Main comparison

1. I downloaded the phylogeny of all birds from `birdtree.org` and took a random sample of 50 phylogenetic trees.
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
2. *Functional traits*: I extracted foraging traits from a database of bird functional traits (EltonTraits, `http://dx.doi.org/10.1890/13-1917.1`). I extracted life-history traits from a database of amniote life history traits (`https://doi.org/10.6084/m9.figshare.c.3308127.v1`). I tested whether the mean difference varied by the literature value of species mean body mass, or by a categorical diet classification. 
3. *Range where the specimens were collected*: I used the coordinates of each point where specimens were collected to calculate two measurements of the spatial spread of the collected specimens. First, I projected the lat-long coordinates into a grid coordinate reference system so that I could get correct areas. I calculated spatial spread as the mean pairwise distance between the coordinates, as well as by calculating the area of the hull containing the points. I tested whether the mean difference was correlated with this spread.
4. *Climate variability where the specimens were collected*: I downloaded the global 1979-2013 BioClim variables for the entire terrestrial globe at 0.5' resolution (`http://dx.doi.org/doi:10.1594/WDCC/CHELSA_v1_1`). I calculated the following variables and tested whether the mean difference was correlated with them: 
	- spatial CV of temperature across the range where specimens were collected
	- spatial CV of precipitation across the range
	- mean seasonal variability of temperature, averaging all the points where specimens were collected (temporal var within year)
	- mean seasonal variability of precipitation (temporal var within year)
	- mean interannual variability of temperature (temporal var across years)
	- mean interannual variability of precipitation (temporal var across years)
5. *Range sizes*: I downloaded the Birds of the World (BotW) range maps from BirdLife International (birdlife.org); I was able to find a range for all species except the green heron which may have been lumped with its sister species--we might have to throw that pair out of the analysis which should not be a problem. I used QGIS to calculate the area of the polygons in square km. I tested whether a species' CV was related to its range size. I used the sum of resident and breeding range as the range size of a species (not counting winter range and areas passed through while migrating).
6. *Migrant status*: I got information on migrant status from the range maps as follows: the range polygons are classified as year-round resident range, breeding range, and nonbreeding (winter) range. I classified species with only a year-round range as nonmigratory, species with both resident and breeding range as partial migrants, and species with no permanent range as obligate migrants. I tested whether the migrant status of a species was related to its CV.
7. *Richness of co-occurring species*: For each geographic coordinate where a specimen was collected, I counted the number of BotW range polygons that overlapped that point to get a rough estimate of richness, or at least the size of the regional species pool. In addition, I counted the number of congeners of the target species whose ranges overlapped the point. Both of these were included as covariates. 

<div style="page-break-after: always;"></div>
## Results

### Key result

Tropical bird species have lower CV of body mass than their temperate sister species: `$t_{100} = -3.06$`, `$p = 0.0015$`, mean difference `$\bar {d} = -0.020$` (less in the tropics), `$CI = (-\infty, -0.0092)$`. (Figure 1: Plot connecting each pair with a line; Figure 2: Histogram of differences).

This result held up to the bootstrapping test, indicating that it is not an artifact of different sample sizes between tropical and temperate birds. The t-test indicated that the tropical species had a significantly lower CV in 95.1% of the bootstrap simulations. The central 95% of the distribution of mean differences from the simulations was `$(-0.024, -0.010)$` and the central 95% of the distribution of t-statistics was `$(-3.22, -1.53)$`.

### Variables that might moderate the relationship

My initial analysis shows that the tropical-temperate difference in variability might be driven by the difference in seasonal (Figure 3) and interannual temperature variability (Figure 4) between the temperate and tropical locations that the specimens were collected. As the difference between temperate and tropical temperature variability increases, the difference between temperate and tropical body mass variability also increases. In addition, the difference is moderated by the average body mass. There is a bigger difference between tropical and temperate in bigger taxa (Figure 5). However, both of these relationships are relatively weak. They only explain a few percent of the variation in each case, with a lot of noise. 

There is also a possible influence of a couple additional factors on this relationship: one is the area from which the specimens were collected. Within species, as specimens were collected from a bigger area, the CV goes up slightly. However, this was essentially random with regard to the temperate-tropical status, so it does not have any relevance for our relationship. And again, it is a weak relationship.

The other covariates (range size, total and congeneric richness, size of the geographic area from which the museum specimens were collected, the spatial variability in climate over which the specimens were collected, the number of populations from which the specimens were collected, and the elevational range from which they were collected) did not show any trend (Figure 6).

The ultimate message here is that we have a net greater amount of variability in temperate taxa, but noise in the relationship that so far we don't have a good way of explaining.  

<div style="page-break-after: always;"></div>
## Figures

![bird figure](file:///C:\\Users\\Q\\Dropbox\\projects\\verts\\vertnet_results\\birdinteractionplot03feb.png)  
**Figure 1.** Interaction plot showing that the tropical species has lower intraspecific variability in body mass than the non-tropical species. Each point is the CV of log-transformed body mass for each species. Sister taxa are connected with a gray line. The red points are the mean CVs for each realm +/- 1 SE, connected by a red line. 

![histogram](file:///C:\\Users\\Q\\Dropbox\\projects\\verts\\vertnet_results\\birdhistogram03feb.png) 
**Figure 2.** Histogram of the mean differences for birds showing they are normally distributed and that the paired t-test is appropriate. The thick red line is the mean difference, the thin red line is the upper bound of the confidence interval, and the dotted blue line is zero. This shows that the mean difference is significantly less than zero.

![seasonal temp scatter](file:///C:\\Users\\Q\\Dropbox\\projects\\verts\\vertnet_results\\temp_seasonal_scatter.png)  
**Figure 3.** Scatter plot of the CV differences plotted against difference in temperature seasonality for each pair. A regression fit is included. Points to the left of the vertical red line have reduced temperature seasonality for the tropical member of the pair (most tropical climates have lower temperature seasonality). Points below the horizontal red line have reduced body size variation for the tropical member of the pair. Janzen's hypothesis would predict a positive trend, which is moderately supported here (`$R^2 = 0.06$`).

![seasonal temp scatter](file:///C:\\Users\\Q\\Dropbox\\projects\\verts\\vertnet_results\\temp_interannual_scatter.png)  
**Figure 4.** Scatter plot of the CV differences plotted against difference in interannual temperature variability for each pair. A regression fit is included. See Figure 3 caption (`$R^2 = 0.13$`).

![seasonal temp scatter](file:///C:\\Users\\Q\\Dropbox\\projects\\verts\\vertnet_results\\log_body_mass_scatter.png)  
**Figure 5.** Scatter plot of the CV differences plotted against mean `$\log_{10} mass$` of each pair. A regression fit is included. Points below the horizontal red line have reduced body size variation for the tropical member of the pair (`$R^2 = 0.07$`).

<img src="file:///C:\\Users\\Q\\Dropbox\\projects\\verts\\vertnet_results\\other_covariates_scatter.png" alt="other covariates scatter" height="662" width="375">  
**Figure 6.** Other proposed covariates that do not seem to influence the within-pair differences.  

