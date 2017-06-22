---
layout: default
---


![bammtrees]({{ site.baseurl }}/assets/whales_sepRateShiftConfigs.png) 

### Overview
Welcome to the 2017 Short Course on Macroevolution organized by [Dan Rabosky](http://www-personal.umich.edu/~drabosky/Home.html) (University of Michigan) and [Brian Sidlauskas](http://people.oregonstate.edu/~sidlausb/) (OSU).  This course will introduce you to several current approaches to inferring macroevolutionary dynamics with phylogenies including [BAMM](http://bamm-project.org/).

#### Location: Room 032 Nash Hall 
[map](http://oregonstate.edu/campusmap/locations/info/821)

#### Course Url: [https://tinyurl.com/osu-macro-2017](https://tinyurl.com/osu-macro-2017)

#### Instructors
[Dan Rabosky](http://www-personal.umich.edu/~drabosky/Home.html), University of Michigan  
[Michael Alfaro](http://pandorasboxfish.squarespace.com/), UCLA  
[Jonathan Mitchell](https://lsa.umich.edu/eeb/people/postdoctoral-fellows/jonsmitc.html), University of Michigan  
[Jonathan Chang](https://jonathanchang.org/), UCLA  
[Pascal Title](http://pascaltitle.weebly.com/), University of Michigan  

### Schedule

#### Thursday 22 June

- R basics. [script]({{ site.baseurl }}/assets/1.general_diversification.R)   *Rabosky*
- Simulating trees with temporal variation in rates. [script]({{ site.baseurl }}/assets/2.temporal_variation.R)  *Rabosky*
- Simulating trees with rate variation among clades. [script]({{ site.baseurl }}/assets/3.among_clade_variation.R)  *Rabosky*
- Trait dependent diversification. [script]({{ site.baseurl }}/assets/4.trait-dependent-diversification.R)  *Rabosky*
- Rate variable models in RPANDA. [link]({{ site.baseurl }}/assets/r-panda.html)  [pdf]({{ site.baseurl }}/assets/r-panda.pdf)  *Alfaro*





#### Friday 23 June 

- fossilBAMM [script]({{ site.baseurl }}fossil-Rscripts/horse_occs.R) *Mitchell*

### Supporting functions

- [diversification functions]({{ site.baseurl }}/assets/1.general_diversification.R)
- [trait dependent functions]({{ site.baseurl }}/supporting/traitDependent_functions.R)
- [shift tree simulator]({{ site.baseurl }}/supporting/simulate_shift_trees.R)

### Datasets 

#### accipitrids

* [accipitrids/accip.csv]({{ site.baseurl }}/data/accipitrids/accip.csv)
* [accipitrids/accipiters.nex]({{ site.baseurl }}/data/accipitrids/accipiters.nex)

#### actinopterygii

* [actinopterygii/bigfish.tre]({{ site.baseurl }}/data/actinopterygii/bigfish.tre)
* [actinopterygii/fish_logsize.txt]({{ site.baseurl }}/data/actinopterygii/fish_logsize.txt)

#### anolis

* [anolis/GA_Anolis_MCC.tre]({{ site.baseurl }}/data/anolis/GA_Anolis_MCC.tre)

#### birds

* [birds/birds_master_taxonomy.csv]({{ site.baseurl }}/data/birds/birds_master_taxonomy.csv)
* [birds/hackett_nozero.tre]({{ site.baseurl }}/data/birds/hackett_nozero.tre)

#### canidae

* [canidae/canidae.tre]({{ site.baseurl }}/data/canidae/canidae.tre)
* [canidae/control.txt]({{ site.baseurl }}/data/canidae/control.txt)
* [canidae/fossil_event_data.txt]({{ site.baseurl }}/data/canidae/fossil_event_data.txt)
* [canidae/fossil_mcmc_out.txt]({{ site.baseurl }}/data/canidae/fossil_mcmc_out.txt)
* [canidae/sampling.txt]({{ site.baseurl }}/data/canidae/sampling.txt)

#### fisse testdata

* [fisse_testdata/example_trait.csv]({{ site.baseurl }}/data/fisse_testdata/example_trait.csv)
* [fisse_testdata/example_tree.tre]({{ site.baseurl }}/data/fisse_testdata/example_tree.tre)

#### horses

* [horses/control.txt]({{ site.baseurl }}/data/horses/control.txt)
* [horses/horse.tre]({{ site.baseurl }}/data/horses/horse.tre)
* [horses/horse_interval.txt]({{ site.baseurl }}/data/horses/horse_interval.txt)
* [horses/interval_times.txt]({{ site.baseurl }}/data/horses/interval_times.txt)

#### jetz birds

* [jetz_etal_post/jetzbirds.tre]({{ site.baseurl }}/data/jetz_etal_post/jetzbirds.tre)
* [jetz_etal_post/jetzbirds_eventdata.txt]({{ site.baseurl }}/data/jetz_etal_post/jetzbirds_eventdata.txt)

#### salamander

* [salamander/control.txt]({{ site.baseurl }}/data/salamander/control.txt)
* [salamander/control_extant.txt]({{ site.baseurl }}/data/salamander/control_extant.txt)
* [salamander/event_data.txt]({{ site.baseurl }}/data/salamander/event_data.txt)
* [salamander/extant.tre]({{ site.baseurl }}/data/salamander/extant.tre)
* [salamander/mcmc_out.txt]({{ site.baseurl }}/data/salamander/mcmc_out.txt)
* [salamander/salamandridae.tre]({{ site.baseurl }}/data/salamander/salamandridae.tre)


#### skinks

* [skinks/skinks216.tre]({{ site.baseurl }}/data/skinks/skinks216.tre)

#### warblers

* [warblers/warbs.tre]({{ site.baseurl }}/data/warblers/warbs.tre)

#### whales

*  [whales/bamm]({{ site.baseurl }}/data/whales/bamm)
*  [whales/divcontrol.txt]({{ site.baseurl }}/data/whales/divcontrol.txt)
*  [whales/myPriors.txt]({{ site.baseurl }}/data/whales/myPriors.txt)
*  [whales/whaletree.tre]({{ site.baseurl }}/data/whales/whaletree.tre)


### Miscellany

- [r reference card]({{ site.baseurl }}/assets/refcard.pdf)
- [r markdown cheat sheet](https://www.rstudio.com/wp-content/uploads/2015/02/rmarkdown-cheatsheet.pdf)
- [dplyr and tidyr reference card](https://www.rstudio.com/wp-content/uploads/2015/02/data-wrangling-cheatsheet.pdf)





