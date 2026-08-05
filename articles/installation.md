# Installation

## Requirements/dependencies

- R (\>= version 4.1)
- **quarto CLI**
  - Downloads available for different platforms here:
    <https://github.com/quarto-dev/quarto-cli/releases>
  - When installing on Linux, make sure you include `quarto` in your
    PATH variable, see instructions
    [here](https://quarto.org/docs/download/tarball.html)
  - Ensure that quarto is correctly installed through `quarto check`

## Installation of oncoEnrichR with R commands

1.  `install.packages('remotes')`
2.  `remotes::install_github('sigven/oncoEnrichR', ref = "v1.6.2")`
3.  [`library(oncoEnrichR)`](https://sigven.github.io/oncoEnrichR)

  
