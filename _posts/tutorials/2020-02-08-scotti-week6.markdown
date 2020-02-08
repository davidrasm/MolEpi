---
layout: page
title: Transmission tree reconstruction with SCOTTI
permalink: /tutorials/scotti-week6/
---

### What you will need:

-[BEAST 2.5.0][beast] or greater <br>
-[Tracer 1.7.0][tracer] or greater <br>
-[FigTree][figtree] <br>
-Python 2.7 or greater

[beast]: <http://www.beast2.org/>
[tracer]: <https://github.com/beast-dev/tracer/releases/tag/v1.7.1>
[figtree]: <http://tree.bio.ed.ac.uk/software/figtree/>

### Transmission tree reconstruction with SCOTTI

In this tutorial we will use SCOTTI (Structured COalescent Transmission Tree Inference) to reconstruct transmission trees for small outbreaks [(De Maio *et al.*, 2016)][demaio-2016]. A great tutorial introducing SCOTTI has already been created by Louis du Plessis and Nicola de Maio, so we will follow along with their tutorial available on the Taming the BEAST website:

* [Link to SCOTTI tutorial][scotti-tutorial]

[demaio-2016]: <https://doi.org/10.1371/journal.pcbi.1005130>
[scotti-tutorial]: <https://taming-the-beast.org/tutorials/SCOTTI-Tutorial/>

I recommend downloading/cloning the entire [repository][git-repo] for the tutorial so you will have all the data files and helper scripts in the same place.

[git-repo]: <https://github.com/taming-the-beast/SCOTTI-Tutorial>


### Reconstructing a neonatal *Klebsiella pneumoniae* outbreak

In addition to the FMDV dataset the Taming that the BEAST tutorial explores, [De Maio *et al.* (2016)][demaio-2016] analyzed an outbreak of antimicrobial resistant *K. pneumoniae* in a Nepali neonatal intensive care unit. This was a severe outbreak that resulted in the deaths of 16 out of 25 infected infants. As an alternative to the FMDV outbreak dataset, you can follow along with the Taming the BEAST tutorial while using the files below to replicate the *K. pneumoniae* analysis. Note, however, that this analysis will take quite a bit longer to run since this it is a larger outbreak with more sampled hosts.

* [Klebsiella sequence data][kleb-seqs] <br>
* [Klebsiella sample dates][kleb-dates] <br>
* [Klebsiella sample hosts][kleb-hosts] <br>
* [Klebsiella host times][kleb-host-times] <br>

[kleb-seqs]: <{{site.baseurl}}/tutorials/scotti-week6/KPneu.fasta>
[kleb-dates]: <{{site.baseurl}}/tutorials/scotti-week6/KPneu_dates.csv>
[kleb-hosts]: <{{site.baseurl}}/tutorials/scotti-week6/KPneu_hosts.csv>
[kleb-host-times]: <{{site.baseurl}}/tutorials/scotti-week6/KPneu_hostTimes.csv>

If you cloned the SCOTTI-Tutorial repo, you can place these files in the *'Data'* folder alongside the FMDV data. You should then be able to run the *SCOTTI_generate_xml.py* as described in the tutorial, but replace the file name arguments with the *K. pneumoniae* file names as in the command below:

```
python SCOTTI_generate_xml.py --fasta ../data/KPneu.fasta --dates ../data/KPneu_dates.csv --hosts ../data/KPneu_hosts.csv --hostTimes ../data/KPneu_hostTimes.csv --output KPneu --maxHosts 40 --numIter 4000000 --tracelog 2000 --treelog 20000 --screenlog 20000
```

Running the python script should generate a *KPneu.xml* file that you can run in BEAST. Note that we've also increased the *--maxHosts* argument to 40 since we know that at least 25 hosts were infected. You should then be able to follow along with the rest of the Taming the BEAST tutorial without any problems.

After processing the BEAST output, I get a MCC tree that looks like this:

<img src="{{site.baseurl}}/assets/img/tutorials/scotti-week6/KPNeu-MCC-hosts.png" alt="Exploring the predictor variables in Tracer" width="700" height="450">

As in [De Maio *et al.* (2016)][demaio-2016], most infections appear to be linked by unsampled hosts.


