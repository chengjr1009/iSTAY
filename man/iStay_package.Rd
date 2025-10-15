\name{iSTAY-package}
\alias{iSTAY-package}
\alias{iSTAY}
\docType{package}
\title{
  Information-based stability and synchrony measures.
}
\description{
  \code{iSTAY} (information-based stability and synchrony measures) is an R package that provides functions to compute a continuum of information-based measures for quantifying the temporal stability of populations, communities, and ecosystems, as well as their associated synchrony, based on species (or species assemblage) biomass or other key variables. When biodiversity data are available, the package also enables the assessment of the corresponding diversity–stability relationships. The information-based measures are derived from Hill numbers parameterized by an order q > 0; see Chao et al. (2025) for the theoretical and methodological background. All measures are illustrated using temporal biomass data from the Jena experiment (Roscher et al. 2004; Weisser et al. 2017; Wagg et al. 2022).  \cr

Specifically, iSTAY provides the  following measures for three types of analyses:

(1) \strong{Single time series}: computes stability measures of order q > 0 and displays the corresponding stability profile for a time series (or for each time series analyzed individually). The stability profile illustrates how stability varies with the order q. When biodiversity data are available, iSTAY also assesses the diversity–stability relationship across individual time series.

(2) \strong{Multiple time series}: computes four measures--gamma, alpha, and beta stability, as well as synchrony--and displays the corresponding profiles for multiple time series (or for each set of time series within a collection). When biodiversity data are available, iSTAY also assesses the diversity–stability and diversity–synchrony relationships.

(3) \strong{Hierarchical Series}: computes four measures--gamma, alpha, and beta stability, as well as synchrony--for each hierarchical level, and provides the corresponding stability and synchrony profiles.\cr


This package contains five main functions: \cr

  1. \code{iSTAY_Single} calculates stability measures of order q > 0 for a single time series of biomass or other relevant variables. \cr

  2. \code{iSTAY_Multiple} computes gamma, alpha, and beta stability, as well as synchrony for multiple time series of biomass or other relevant variables. \cr

  3. \code{iSTAY_Hier} computes gamma, alpha, and beta stability, as well as synchrony, at each hierarchical level for hierarchical time series of biomass or other variables. \cr

  4. \code{ggiSTAY_qprofile} plots the stability and synchrony profiles based on the output obtained from the functions \code{iSTAY_Single}, \code{iSTAY_Multiple} or \code{iSTAY_Hier}. \cr

  5. \code{ggiSTAY_analysis} plots the diversity–stability and diversity–synchrony relationships based on the output obtained from the function \code{iSTAY_Single} or \code{iSTAY_Multiple}. \cr

}

\references{
Chao, A., Colwell, R. K., Shia J., Thorn, S., Yang, M.-Y., Mitesser, O., et al. (2025) A continuum of information-based temporal stability measures and their decomposition across hierarchical level. BioRxiv < https://doi.org/10.1101/2025.08.20.671203>\cr\cr

Roscher C. Schumacher, J., Baade, J., Wilcke, W., Gleixner, G., Weisser, W. W. et al. (2004). The role of biodiversity for element cycling and trophic interactions: an experimental approach in a grassland community. Basic and Applied Ecology, 5, 107–121.\cr\cr

Wagg, C., Roscher, C., Weigelt, A., Vogel, A., Ebeling, A., De Luca, E. et al. (2022) Biodiversity–stability relationships strengthen over time in a long-term grassland experiment. Nature Communications, 13, 7752.\cr\cr

Weisser, W. W., Roscher, C., Meyer, S. T., Ebeling, A., Luo, G., Allan, E. et al. (2017). Biodiversity effects on ecosystem functioning in a 15-year grassland experiment: Patterns, mechanisms, and open questions. Basic and Applied Ecology, 23, 1–73.\cr\cr

}

\keyword{ package }
