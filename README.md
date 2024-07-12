# survivl

This package allows for simulation of survival data from additive hazard and Cox Marginal Structural Models.

Current features:

-   Simulation of Cox MSM survival models.

To install and load the package, run the commands

```         
install.packages("devtools")
devtools::install_github("rje42/causl") # dependency
devtools::install_github("rje42/survivl")
library(survivl)
```

## Example:

```mermaid
graph TD;
    A[Z_{L-1}] --> B[X_{L-1}];
    B --> C[Y_{L-1}];
    A --> C;
    B --> D;
    D[Z_{L}] --> E[X_{L}];
    E --> F[Y_{L}];
    D --> F;
