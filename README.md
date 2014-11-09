qEfficientFrontier
==================

R package for QlikView that takes in any portfolio of stocks return and optimizes it using the Efficient Frontier
This is the r package used with the qEfficientFrontier to process data send from a QlikView Document to R and back. The R package contains only one function called eff_frontier which takes a json object called Data. The R package is build in a way so that it will always send the exact data structure back to QlikView with the new allocation weights.
