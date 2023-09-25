Use CDO grid files and template shell for generating grist accepted initial data.
The initialization procedure depends on de vertical level and coordinate, besides
a horizontal remapping.

Current code support:
1. ERAIM-60: ERA-interim model-level data   , 60 levels
1. ERAIP-37: ERA-interim pressure-level data, 37 levels

For other data source, the current procedure needs to be consolidated that
it can be entirely generalzied, an issue that has not been examined.
ERA-interim and ERA-5 are recommended.
