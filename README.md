# maicChecks <a href="https://github.com/clyau/maicChecks"><img src="https://github.com/user-attachments/assets/16780205-98c6-4bbd-af23-dfec94010d54" align="right" height="140"/></a>

<!-- badges: start -->

[![CRAN status](https://badges.cranchecks.info/flavor/release/maicChecks.svg)](https://cran.r-project.org/web/checks/check_results_maicChecks.html) [![downloads](https://cranlogs.r-pkg.org/badges/maicChecks)](https://www.rdocumentation.org/trends)

<!-- badges: end -->

`maicChecks`: This package provides tools for data analysts to conduct additional statistical checks and ensure the numerical feasibility of implementing MAIC, a method widely used in clinical trials and observational studies for comparing treatment effects.

## Installation

| Type                | Source | Command                                       |
|------------------|------------------|------------------------------------|
| Release (0.1.2)     | CRAN   | `install.packages("maicChecks")`              |
| Development (0.2.0) | GitHub | `remotes::install_github("clyau/maicChecks")` |

## Overview

The Geometric Approaches to Assessing Numerical Feasibility for Conducting MAIC package offers a comprehensive set of tools to assess the appropriateness of MAIC for a given set of baseline characteristics. The package includes methods for checking the overlap between individual patient data (IPD) and aggregated data (AD), visualizing data distributions, and performing statistical tests to ensure valid comparisons.

## Key Features

-   **Convex Hull Check**: Ensures that the AD lies within the convex hull of the IPD, guaranteeing a unique solution for MAIC weights.
-   **Principal Component Analysis (PCA)**: Provides a visual assessment of the AD's position relative to the IPD in a multi-dimensional space.
-   **Mahalanobis Distance and Hotelling's T² Test**: Tests whether matching IPD to AD is necessary by assessing the similarity of their distributions.
-   **Exact Matching**: An alternative to propensity score matching, offering exact matching of covariate means between treatment groups.

## Methods

-   **Convex Hull Check**: This method uses linear programming to determine if the AD is within the convex hull of the IPD, ensuring numerical compatibility for MAIC.
-   **Principal Component Analysis (PCA)**: PCA is used to visualize the AD's position relative to the IPD, providing a graphical representation of data overlap.
-   **Mahalanobis Distance and Hotelling's T² Test**: These statistical tests assess the similarity between IPD and AD, determining if matching is necessary.
-   **Exact Matching**: This method provides exact matching of covariate means between treatment groups, ensuring balanced comparisons.
