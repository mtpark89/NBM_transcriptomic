# NBM expression project
Code to extract and characterize brain region specific expression (primarily from the Allen Brain Institute atlases)

This code corresponds to the article titled: "Transcriptomic Characterization of the Human Insular Cortex and Claustrum"
https://www.frontiersin.org/articles/10.3389/fnana.2019.00094/full

This is a rough fork from the https://github.com/leonfrench/HumanHabenula project that generalized it and applied it to the basal forebrain regions (diagonal band and nucleus basalis of Meynert).


### Steps
* download the Allen Adult human brain atlas to /data/raw/allen_HBA
* download the Allen+Brainspan fetal human brain atlas to /data/raw/allen_human_fetal_brain
* run the R scripts in numbered order

### References for data sources used

Mancarci, Ogan B., Lilah Toker, Shreejoy J. Tripathy, Brenna Li, Brad Rocco, Etienne Sibille, and Paul Pavlidis. 2017. “Cross-Laboratory Analysis of Brain Cell Type Transcriptomes with Applications to Interpretation of Bulk Tissue Data.” eNeuro, November, ENEURO.0212–17.2017.

Miller, Jeremy A., Song-Lin Ding, Susan M. Sunkin, Kimberly A. Smith, Lydia Ng, Aaron Szafer, Amanda Ebbert, et al. 2014. “Transcriptional Landscape of the Prenatal Human Brain.” Nature 508 (7495): 199–206.

Hawrylycz, M. J., Lein, E. S., Guillozet-Bongaarts, A. L., Shen, E. H., Ng, L., Miller, J. A., et al. (2012). An anatomically comprehensive atlas of the adult human brain transcriptome. Nature 489, 391–399.

Darmanis, S., Sloan, S. A., Zhang, Y., Enge, M., Caneda, C., Shuer, L. M., et al. (2015). A survey of human brain transcriptome diversity at the single cell level. Proc. Natl. Acad. Sci. U. S. A. 112, 7285–7290.

Piñero, Janet, Àlex Bravo, Núria Queralt-Rosinach, Alba Gutiérrez-Sacristán, Jordi Deu-Pons, Emilio Centeno, Javier García-García, Ferran Sanz, and Laura I. Furlong. 2017. “DisGeNET: A Comprehensive Platform Integrating Information on Human Disease-Associated Genes and Variants.” Nucleic Acids Research 45 (D1): D833–39.
