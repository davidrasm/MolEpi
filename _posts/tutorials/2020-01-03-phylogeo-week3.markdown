---
layout: page
title: Discrete trait phylogeography
permalink: /tutorials/phylogeo-week3/
---

### What you will need:

-[BEAST v.2.5.0][beast] or greater <br>
-[Tracer v1.7.0][tracer] or greater <br>
-[FigTree][figtree] <br>

[beast]: <http://www.beast2.org/>
[tracer]: <https://github.com/beast-dev/tracer/releases/tag/v1.7.1>
[figtree]: <http://tree.bio.ed.ac.uk/software/figtree/>

### Exploring the origins of *Phytophthora infestans*

In this tutorial, we will use discrete trait phylogeographic methods to explore the origins *Phytophthora infestans*, the pathogen behind the Irish potato famine. Two origins have been postulated for the origins of P. infestans: one out of the Andean region of South America and one out of central Mexico. 

The data for this tutorial was kindly provided by Niklaus Grunwald and Javier Tabima at Oregon State.

### Phylogeography in BEAST 2

The discrete trait phylogeographic models we will use in this tutorial were originally implemented in BEAST v1, but can be accessed in BEAST 2 through the ***BEASTlabs*** and ***BEAST-CLASSIC*** add on packages. 

To install these packages, open BEAUTi and select ***File &rarr; Manage Packages***. In the dialog box that pops up, first select BEASTlabs and then click the Install button. Then install the BEAST-CLASSIC packages. Once both packages are installed, ***restart BEAUTi*** so that the changes take affect in BEAUTi. 

<img src="{{site.baseurl}}/assets/img/tutorials/week3/install_BEASTlabs.png" alt="Installing BEASTlabs and BEAST-CLASSIC" width="600" height="400">

Reopen BEAUTi and click the ***+*** button on bottom right of the Partitions tab. Under select what to add select ***Import Alignment*** and click OK. Then import the ***phyt_IRRAS.fasta*** aligment file. 

<img src="{{site.baseurl}}/assets/img/tutorials/week3/import_align.png" alt="Importing the alignment in BEAUTi" width="600" height="400">

<img src="{{site.baseurl}}/assets/img/tutorials/week3/imported_phyt_align.png" alt="The imported alignment" width="600" height="400">

We will assume all samples are contemporaneous (i.e. sampled at the present) so we do not need to add anything in the Tip Dates panel.

In the Site Model panel, we will use a gamma site model with four rate categories, so enter '4' next to Gamma Category Count. A HKY substitution model should be sufficient for this data. The Site Model panel should look like the screenshot below.

<img src="{{site.baseurl}}/assets/img/tutorials/week3/site_model.png" alt="Setting up the site model" width="600" height="400">

In the Clock Model panel we do not need to change anything because we will use a strict clock model and we will not estimate the Clock.rate since we do not have a way of calibrating the clock rate for the RAS gene. Time will therefore be in units of arbitrary time for our analysis, but we can still reconstruct the ancestral locations of lineages through time. BEAST 2 however likes to force you to either estimate the clock rate or the substitution rate, so to prevent this from happening click ***Mode*** and then unselect ***Automatic set clock rate***.

<img src="{{site.baseurl}}/assets/img/tutorials/week3/clock_model.png" alt="Unselecting the automatic clock rate setting" width="600" height="400">

In the Priors panel, we will use a ***Coalescent Constant Population*** tree prior. The other priors can be left at their default settings.

<img src="{{site.baseurl}}/assets/img/tutorials/week3/tree_prior.png" alt="Settting the coalescent tree prior" width="600" height="400">

Now we want to add our data about our discrete trait, which will be the geographic location of each sample. BEAST 2 treats this trait data in a very similar way to adding another data set or partition to the alignment. Go back to the Partitions screen. 

---
**Note**

In the orignal analysis of Goss et al. (2014), 13 discrete states were used. While it is computationally feasible to have this many locations, this means that there will be 13 x 13 = 169 migration rate parameters need to be estimated.  So to simplify things, I've grouped all countries in continental Europe (Estonia,Hungary,Netherlands,Poland,Sweden) into one location 'EU'. Furthermore I eliminated two countries (Vietnam and South Africa). This results in a 7 x 7 rate matrix. Grouping or binning nearby locations is a common strategy.

---

### Acknoweledgments

This tutorial was based on an [earlier tutorial][remco_tutorial] by Remco Bouckaert on discrete traits models in BEAST 2.

[remco_tutorial]: <https://github.com/BEAST2-Dev/beast-classic/releases/download/v1.3.0/ARv2.4.1.pdf>






