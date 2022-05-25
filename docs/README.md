# XRF-Maps

## About

## Fitting Routines
```
--fit <roi,nnls,svd,matrix>
```
### ROI
#### Description
ROI fitting 

### Non-Negative Least Squares (NNLS)

### Singualr Value Decomposition (SVD)

### Per Pixel (MATRIX)

## Quantification
### Option
```
--quantify-with <filename.txt>
```
### Standard file definition:

### Example:
```
--quantify-with maps_standardinfo.txt
```

## Exchange Format

## MISC Command line 

### Number of Threads
```
--nthreads <int>
```
#### Desc: 
Specifiy the max number of threads to use. If you do not specifiy --nthreads, it will query the cpu for the max number of threads and use them all.

#### Example: Use only 4 threads
```
--nthreads 4
```

### Detectors to fit
```
--detectors <n,...,m>
```
or
```
--detector-range n:m
```
#### Desc: 
Specifiy the detector elements to fit. If you do not specifiy this option, it will default to 7 elements but will only process ones that exist.

#### Example: Use 4 detectors
```
--detectors 0,1,2,3
```
or
```
--detector-range 0:3
```
#### Example: Fit all detectors except for the 3rd detector
```
--detectors 0,1,3
```


### Generate Average H5 file
```
--generate-avg-h5
```
#### Desc: 
Will average all detectors specified by --detectors or --detector-range option into h5 file. 





<!-- ![This is an image](https://myoctocat.com/assets/images/base-octocat.svg) -->