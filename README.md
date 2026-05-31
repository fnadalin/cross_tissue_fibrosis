# Cross-tissue cell neighbourhood analysis

This repository collects code for comparing disease-associated cell neighbourhoods across tissues.
The approach consists of two main steps:

* Per-tissue step:
  * computation of disease differentially abundance neighbourhoods with Milo;
  * neighbourhood aggregation into neighbourhood groups (nhg);
  * calculation of the abundance logFC between disease and normal and the markers of nhgs.
* Cross-tissue step:
  * construction of a multi-partite graph G where each tissue is a layer and each nhg is a node;
  * define a weighted version G' of G, where each edge weight is a function of the abundance logFC and of the Jaccard index of the markers of the adjacent nodes;
  * clustering of the multi-partite weighted graph G'.

