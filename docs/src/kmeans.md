# K-Means Clustering

K-means clustering of multidimensional vector data based on L2 distances is implemented. No specialised data structures are utilised and data is expected to be stored in matrices of Floats. The methods also do not check dimensions of inputs or bounds in for loops - *it is expected that users will take care to avoid mistakes*. Currently the emphasis is on speed (not convenience) with a view to implement GPU versions of these methods in the future.

## Clustering Methods

```@docs
kmeans
kmeans!
```


## Support Methods
The following methods that support the kmeans algorithm implementation are also exposed to users

```@docs
kmeans_init
closest_centre
closest_centre!
mindist2cntrs
```