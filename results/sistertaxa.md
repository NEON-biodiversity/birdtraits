# Variation in sister species

QDR, 19 Jan 2017

## Summary

The prediction here is based off the general dogma that species have narrower ranges/niches in the tropics relative to places outside the tropics, because they are adapted to the narrower range of environmental variability in the tropics. That leads us to predict that for a pair of closely related species in which one is in the tropics and one is not, the tropical species will have less trait variability. I tested this prediction using body mass measurements on vertebrate museum specimens from vertnet.org. The prediction is supported for birds! It is not supported for mammals but there are many fewer pairs to compare. I could not find enough data to test the prediction in amphibians, fish, or reptiles.

## Methods

1. For each group of vertebrates, pull all records from vertnet.org that have body mass measured.
2. Find genera that have multiple species with at least 5 body mass records.
3. Find the median latitude of all the records for each species.
4. Find pairs within genera in which the lower-latitude species is tropical (`$|lat| < 23.5$`) and the higher-latitude species is not tropical (`$|lat| > 23.5)$`.  
	- For birds, there are 98 such pairs with adequate number of records.
	- For mammals, there are 15 pairs.
	- For reptiles, there are only 2, and there are 0 for amphibians and fish, so the remaining analysis will focus on birds and mammals.
5. Calculate the coefficient of variation of `$\log_{10} mass$` for each species. This is the ITV value for that species.
6. Run a paired one-sided t-test (the most powerful statistical test currently existing today!) using the alternative hypothesis that the tropical species will have lower CV.

## Results

The prediction is supported in birds! The tropical species have lower CV of body mass. (If the taxa are all pooled and a single t-test is run, it is also significant, but since birds are 98/115 of the pairs, that is just reproducing the bird result.) In mammals, the mean difference between temperate and tropical CV was not different from 0.

**Birds:** `$t_{97} = -4.41$`, `$p = 1.4 \times 10^{-5}$`, mean difference -0.049 (less in the tropics), CI `$(-\infty, -0.031)$`  
**Mammals:** `$t_{14} = 0.35$`, `$p = 0.63$`, mean difference 0.022 (greater in the tropics), CI `$(-\infty, 0.13)$`

![bird figure](file:///C:\\Users\\Q\\Dropbox\\neon\\data\\vertnet_results\\sistertaxa_birds.png)  
**Figure 1.** Body mass variability of bird sister species pairs in and out of the tropics. The mean and standard error of each group is in red.

![mammal figure](file:///C:\\Users\\Q\\Dropbox\\neon\\data\\vertnet_results\\sistertaxa_mammals.png)  
**Figure 2.** Body mass variability of mammal sister species pairs in and out of the tropics. The mean and standard error of each group is in red.

![mammal figure](file:///C:\\Users\\Q\\Dropbox\\neon\\data\\vertnet_results\\hist_meandiffs.png)  
**Figure 3.** Histogram of the mean differences for birds showing they are normally distributed and that the paired t-test is appropriate. The thick red line is the mean difference, the thin red line is the upper bound of the confidence interval, and the dotted blue line is zero. This shows that the mean difference is significantly less than zero.