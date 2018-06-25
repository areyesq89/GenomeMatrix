# GenomeMatrix
In the analysis of HiC or similar types of interaction data, somteimes it
is useful to plot matrices rotating them by -45 degrees. My motivation to
write this function was to be able to visualize HiC data along with other
genomic tracks, such as gene annotations (using ggbio) and chromatin marks.
I wanted to plot this using ggplot2, given all its fancy features and
advantages. I show an example below of the resulting plots:

![alt text](https://pbs.twimg.com/media/DgVB7cMWsAE8Eec.jpg)


## Installing the package

```
devtools::install_github("areyesq89/GenomeMatrix")
```