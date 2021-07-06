# dbrpart

Performs CART with distance-based partitions

In its present form, the `dbrpart` package is a prototype that, while working, requires substantial cleanup. This will eventually be submitted to CRAN once this process is complete.

# Missing features

* At present, only supports distance-based splits. This should eventually also support splits on categorical and scalar data
* At present, only supports categorical and distance-based targets. This should eventually also support splits on numeric targets.


# Known Issues

* Optimising combinations of values is implemented using Self-Organising Migrating Algortithm (SOMA). The current implementation and default settings for SOMA make it too slow to practically run, and it might be that some other optimisation algorithm is better suited. More work is needed here.
  * The present workaround is to define combinations of features with different weights and perform a grid search instead
