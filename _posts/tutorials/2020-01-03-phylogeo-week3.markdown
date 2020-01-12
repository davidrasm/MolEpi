---
layout: page
title: Discrete trait phylogeography
permalink: /tutorials/phylogeo-week3/
---

### What you will need:

-[BEAST v.2.5.0][beast] or greater <br>
-[Tracer v1.7.0][tracer] or greater <br>
-[FigTree][figtree] <br>

For postprocessing (optional):
-[Spread][spread] <br>
-[Google Earth Pro][gearth]

[beast]: <http://www.beast2.org/>
[tracer]: <https://github.com/beast-dev/tracer/releases/tag/v1.7.1>
[figtree]: <http://tree.bio.ed.ac.uk/software/figtree/>
[spread]: <https://rega.kuleuven.be/cev/ecv/software/spread>
[gearth]: <https://www.google.com/earth/versions/#earth-pro>

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

Note may want to give a more reasonable rate like 0.00001

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

We can summarize the posterior distribution of trees and reconstructed ancestral states as a single concensus (MCC) tree in TreeAnnotator. Open TreeAnnotator and for the Input Tree File select the *location_tree_with_trait.trees* file that BEAST generated. For the Burnin percentage, we can use 10%. Make sure that Maximum clade credibility tree is selected for Target tree type. Select an Output File name for the MCC tree such as ***phyt_PITG_mcc.tre***. Now press ***Run***.

<img src="{{site.baseurl}}/assets/img/tutorials/week3/summarizing_trees.png" alt="Creating a MCC consensus tree in TreeAnnotator" width="600" height="400">  

Open the consensus MCC tree you just generated in in FigTree. Select ***Appearance***, and then choose location under Colour by. The lineages in the tree will now be colored by their most probable ancestral location (the state with the highest posterior probability). You may want to increase the Line Width to make this easier to see. You can also have the branch widths reflect the posterior support for the ancestral states by selecting location.prob under Width by. 

<img src="{{site.baseurl}}/assets/img/tutorials/week3/mcc_tree_colored.png" alt="The MCMC coloured by ancestral locations" width="600" height="400">

### Visualizing geographic spread

***Note:*** Update Spread link to 1.07.

Download Spread if you haven't already. Unfortunately, the Mac version does not seem to work on newer OS versions, so download the .jar file. Launch it through the command line:

```
java -jar SPREAD\ v1.0.7.jar
```

In Spread, select ***Open &rarr; Load tree file*** and then choose the consensus MCC tree file. Under State attribute name select *location*. 

If you click the ***Set-up*** button, you can manually edit the longitude and latitude locations of each discrete sampling location. I used the cooridinates of the centroid of each country available [here][worldmap-centroids]. 

[worldmap-centroids]: <https://worldmap.harvard.edu/data/geonode:country_centroids_az8>

| Location | Code | Latitude | Longitude |
| ---- | ---- | ---- | ---- |
| Europe | EU | 51.1069 | 10.385 |
| Peru | PE | -9.152 | -74.382 |
| United Kingdom | UK | 54.123 | -2.865 |
| Mexico | MX | 23.947 | -102.523 |
| Columbia | C0 | 3.913 | -73.081 |
| Ecuador | EC | -1.423 | -78.752 |
| United States | US | 45.679 | -112.461 |

---

***Note:*** If you need to find the geographic coordinates of more specific locations, you can always use Google Earth. Pointing the mouse cursor at a specifc location will provide the longitude and latitude.

___

You can then click the output button to generate a KML file. Clicking display will show the migration events in the MCC phylogeny mapped onto the world map. The diameter of the circles shows how many lineages in the tree reside in a particular location at a particular time. The branch lengths are colored by their timing. 

<img src="{{site.baseurl}}/assets/img/tutorials/week3/phyt_global_spread.png" alt="Visualizing the spread of P. infestans in Spread" width="600" height="400">

Finally, you can create an animated video of the spread of lineages through time and space that can be viewed in Google Earth. Download [Google Earth Pro][gearth] and open the application. Click ***File &rarr; Open***, select the KML file you just generated and then click ***Open***. The animation should play automatically and you can rewind or advance through the animation using the cursor at the top of the window.


### Other phylogeography tutorials

Phylogeographic methods and related visualization tools have developed incredibly rapidly in the past decade. So in some sense this tutorial is already out-of-date relative to the cutting edge. So below I've provided links to other tutorials using newer methods. Note that these tutorials are based on BEAST v1.

There is another great tutorial on discrete-trait phylogeography exploring the spread of bat rabies [here][bat-tutorial]. This tutorial is also really neat because it shows how to create web-based visualizations in D3 using SpreaD3, the succesor to Spread, and identify predictors of migration rates using Generalized Linear Models.

[bat-tutorial]: <https://beast.community/workshop_discrete_diffusion>

There is also a very cool tutorial on tracking the spread of West Nile Virus through North America using continuous-trait phylogeographic models [here][continuous-trait-tutorial]. 

[continuous-trait-tutorial]: <https://beast.community/workshop_continuous_diffusion_wnv>

---
**Note**

In the orignal analysis of Goss et al. (2014), 13 discrete states were used. While it is computationally feasible to have this many locations, this means that there will be 13 x 13 = 169 migration rate parameters need to be estimated.  So to simplify things, I've grouped all countries in continental Europe (Estonia,Hungary,Netherlands,Poland,Sweden) into one location 'EU'. Furthermore I eliminated two countries (Vietnam and South Africa). This results in a 7 x 7 rate matrix. Grouping or binning nearby locations is a common strategy.

---

### Acknoweledgments

This tutorial was based on an [earlier tutorial][remco_tutorial] by Remco Bouckaert on discrete traits models in BEAST 2.

[remco_tutorial]: <https://github.com/BEAST2-Dev/beast-classic/releases/download/v1.3.0/ARv2.4.1.pdf>






