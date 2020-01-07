---
layout: page
title: Discrete trait phylogeography
permalink: /tutorials/phylogeo-week3/
---

### What you will need:

-[BEAST v.2.5.0][beast] or greater <br>
-[Tracer v1.7.0][tracer] or greater <br>
-[FigTree][figtree] <br>
-[Spread][spread] <br>

[beast]: <http://www.beast2.org/>
[tracer]: <https://github.com/beast-dev/tracer/releases/tag/v1.7.1>
[figtree]: <http://tree.bio.ed.ac.uk/software/figtree/>
[spread]: <https://www.kuleuven.be/aidslab/phylogeography/SPREAD.html>

### Exploring the origins of *Phytophthora infestans*

In this tutorial, we will use discrete trait phylogeographic methods to explore the origins *Phytophthora infestans*, the pathogen behind the Irish potato famine. Two origins have been postulated for the origins of P. infestans: one out of the Andean region of South America and one out of central Mexico. 

The data for this tutorial was kindly provided by Niklaus Grunwald and Javier Tabima at Oregon State.

### Phylogeography in BEAST 2

The discrete trait phylogeographic models we will use in this tutorial were originally implemented in BEAST v1, but can be accessed in BEAST 2 through the ***BEASTlabs*** and ***BEAST-CLASSIC*** add on packages. 

To install these packages, open BEAUTi and select ***File &rarr; Manage Packages***. In the dialog box that pops up, first select BEASTlabs and then click the Install button. Then install the BEAST-CLASSIC packages. Once both packages are installed, ***restart BEAUTi*** so that the changes take affect in BEAUTi. 

<img src="{{site.baseurl}}/assets/img/tutorials/week3/install_BEASTlabs.png" alt="Installing BEASTlabs and BEAST-CLASSIC" width="600" height="400">

### Setting up the analysis

Reopen BEAUTi and click the ***+*** button on bottom left of the Partitions tab. Select ***Import Alignment***, click OK, and then import the ***phyt_IRRAS.fasta*** alignment file or one of the other four nuclear genes. 

<img src="{{site.baseurl}}/assets/img/tutorials/week3/import_align.png" alt="Importing the alignment in BEAUTi" width="600" height="400">

<img src="{{site.baseurl}}/assets/img/tutorials/week3/imported_phyt_align.png" alt="The imported alignment" width="600" height="400">

We will assume all samples are contemporaneous (i.e. sampled at the present) so we do not need to add anything in the Tip Dates panel.

In the Site Model panel, we will use a gamma site model with four rate categories, so enter '4' next to Gamma Category Count. A HKY substitution model should be sufficient for this data. The Site Model panel should look like the screenshot below.

<img src="{{site.baseurl}}/assets/img/tutorials/week3/site_model.png" alt="Setting up the site model" width="600" height="400">

In the Clock Model panel we do not need to change anything because we will use a strict clock model but we will not estimate the Clock.rate since we do not have a way of calibrating the molecular clock for the nuclear genes used here. Time will therefore be in units of arbitrary time for our analysis, but we can still reconstruct the ancestral locations of lineages through. BEAST however likes to force you to either estimate the clock rate or the substitution rate, so to prevent this from happening click ***Mode*** and then unselect ***Automatic set clock rate***. In the Create new trait windown you can give the trait a name like *location*, then select OK.

<img src="{{site.baseurl}}/assets/img/tutorials/week3/clock_model.png" alt="Unselecting the automatic clock rate setting" width="600" height="400">

In the Priors panel, we will use a ***Coalescent Constant Population*** tree prior. The other priors can be left at their default settings.

<img src="{{site.baseurl}}/assets/img/tutorials/week3/tree_prior.png" alt="Settting the coalescent tree prior" width="600" height="400">

Now we want to add our data about the sampling location of each tip, which will be the discrete trait in our analysis. BEAST treats this trait data very similar sequence data -- in essence we are adding another parition to our sequence alignment with one "site" representing the sampling location of each tip. Go back to the Partitions panel and click the ***+*** button on bottom left. This time select ***Add Discrete Trait***. You can now give the trait a name such as *location*, and then associate the tree with the tree from the sequence alignment.

In the new window that pops up we will now associate each taxon name with its sampling location. Lucky for us, the sampling locations are already contained within the names. Click ***Guess*** and then select use everyting after first underscore. The traits should now be automatically populated with the correct sampling locations, but make sure that the locations correspond to the correct sampling locations in the taxon names!

<img src="{{site.baseurl}}/assets/img/tutorials/week3/importing_traits.png" alt="Importing the discrete trait data" width="600" height="400">

Go to the Site Model panel and select location from the different partitions. The site model here is analagous to the substitution model we used for the sequence data, except now the substitution rate matrix gives the relative migration rate between each location. We will use the default Symmetric model.

<img src="{{site.baseurl}}/assets/img/tutorials/week3/trait_site_model.png" alt="Setting up the trait site model" width="600" height="400">

Go to the Clock Model panel and again select the location partition. We will again use a strict clock model, which controls the absolute rate at which lineages transition between different locations. Make sure the estimate box is checked on the right.

<img src="{{site.baseurl}}/assets/img/tutorials/week3/trait_clock_model.png" alt="Setting up the trait site model" width="600" height="400">

Now we need to set up some priors. For the nonZeroRates we can continue to use a Poisson prior with lambda = 0.693. This prior is used for Bayesian Stochastic Search Variable Selection. Lambda is the prior probability that a given migration rate between two locations will be nonzero. This means that some of the migration rate priors will be constrained to zero, so that the model simplifies to a lower dimensional model with fewer parameters to estimate.

<img src="{{site.baseurl}}/assets/img/tutorials/week3/prior_nonZeroRate.png" alt="Setting the nonZeroRate prior" width="600" height="400">

For the relativeGeoRates and traitClockRate, I recommend Log Normal priors. When in doubt, a Log Normal prior is good prior for rate parameters because the distribution is constrained to positive values. Set the mean M equal to 1.0 and the standard deviation S to 2.0, so that we have a fairly uninformative prior.

<img src="{{site.baseurl}}/assets/img/tutorials/week3/prior_relativeGeoRates.png" alt="Setting the relativeGeoRates prior" width="600" height="400">

Go to the MCMC panel. For a definitive analysis, we would probably want to run the MCMC for much longer, but we don't have all day, so set the Chain Length to 3 million. I would also recommend going under the ***treelog*** drop down menu and setting Log Every to 10,000, otherwise BEAST will write out the tree every 1,000 iterations and the tree file will become huge.

<img src="{{site.baseurl}}/assets/img/tutorials/week3/mcmc_settings.png" alt="Setting the relativeGeoRates prior" width="600" height="400">  

Finally, select ***File &rarr; Save As*** to generate the XML file to input into BEAST. We've now set up our entire analysis! Open BEAST and choose the XML file we just generated and hit ***Run***.

### Analyzing the results

Open Tracer and then import the .log file from your BEAST run. We can use Tracer to explore the posterior estimates of some key parameters like the relGeoRates and the traitClockRate.

### Exploring the phylogeographic reconstructions

We can summarize the posterior distribution of trees and ancestral locations as a single concensus (MCC) tree in TreeAnnotator. Open TreeAnnotator and select the *location_tree_with_trait.trees* file that BEAST created for the Input Tree File. For the Burnin percentage, we can use 10%. Make sure that Maximum clade credibility tree is selected for Target tree type. Select an Output File name for the MCC tree such as ***phyt_PITG_mcc.tre***.  Now press ***Run***.

<img src="{{site.baseurl}}/assets/img/tutorials/week3/summarizing_trees.png" alt="Creating a MCC consensus tree in TreeAnnotator" width="600" height="400">  

Open the consensus MCC tree you just generated in in FigTree. Select ***Appearance***, and then choose location under Colour by. The lineages in the tree will now be colored by their most probable ancestral location (the state with the highest posterior probability). You may want to increase the Line Width to make this easier to see. You can also have the branch widths reflect the posterior support for the most probable ancestral state by selecting location.prob under Width by. 

<img src="{{site.baseurl}}/assets/img/tutorials/week3/mcc_tree_colored.png" alt="The MCMC coloured by ancestral locations" width="600" height="400">

### Visualizing geographic spread

Download Spread if you have not done so already. In Spread, select ***Open &rarr; Load tree file*** and then choose the consensus tree file we generated.

If you click the ***Set-up*** button, you can manually edit the latitude and longitude. locations 

To do: Insert table

### Other tutorials

Continuous phylogeography.

---
**Note**

In the orignal analysis of Goss et al. (2014), 13 discrete states were used. While it is computationally feasible to have this many locations, this means that there will be 13 x 13 = 169 migration rate parameters need to be estimated.  So to simplify things, I've grouped all countries in continental Europe (Estonia,Hungary,Netherlands,Poland,Sweden) into one location 'EU'. Furthermore I eliminated two countries (Vietnam and South Africa). This results in a 7 x 7 rate matrix. Grouping or binning nearby locations is a common strategy.

---

### Acknoweledgments

This tutorial was based on an [earlier tutorial][remco_tutorial] by Remco Bouckaert on discrete traits models in BEAST 2.

[remco_tutorial]: <https://github.com/BEAST2-Dev/beast-classic/releases/download/v1.3.0/ARv2.4.1.pdf>






