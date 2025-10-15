\name{Data_Jena_hierarchical_structure}
\alias{Data_Jena_hierarchical_structure}
\docType{data}
\title{
Four-level hierarchical structure of the Jena Experiment
}
\description{
  All plots in the Jena Experiment were arranged into 4 blocks, each containing 18–20 plots. Plots were sown along a gradient of 1, 2, 4, 8, or 16 species, comprising a total of 462 populations. The biomass of each species sown in a plot was recorded annually from 2003 to 2024, except in 2004. The entire dataset forms a four-level hierarchy: Level 1: population, Level 2: community, Level 3: block, and level 4: overall data. 
This dataset describes the four-level structure. Each row corresponds to a population identifier ( matching the same row \code{Jena_hierarchical_data}), and records the names of block, plot and species to which that population belongs.
}
\usage{data(Data_Jena_hierarchical_structure)}
\format{
  \code{Data_Jena_hierarchical_structure} is a population-by-level (462 × 3) data frame with population identifiers as rows and hierarchical level ("block", "plot" and "species") as columns.
 }

\keyword{datasets}
