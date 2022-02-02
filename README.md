# HT-CETSA-aRPPA-Analysis

The following Matlab script file is provided as open source software under the conditions set forth by 'The MIT License' (see attached license file for details).

```
a190606_PlateSpot_v08_v02b_ForPublication220202.m
```

## Brief description of file function
This matlab script runs in Matlab 2021a and requires the 'image analysis toobox'. The script functions as described in Material and Methods section of peer reviewed scientific journal article title: 
High-throughput cellular thermal shift assay (CETSA) using acoustic transfer of protein lysates
Authors:
Ashley E. Owens1, Michael J. Iannotti1, Tino W. Sanchez1, Ty Voss1, Abhijeet Kapoor1, Matthew D. Hall1, Juan J. Marugan1, Sam Michael1, Noel Southall1, Mark J. Henderson1* 
Affiliation:
1 National Center for Advancing Translational Sciences, National Institutes of Health, Rockville, Maryland, 20850, USA
*Corresponding Author:  Mark Henderson; email: mark.henderson2@nih.gov

Published in ACS Chemical Biology 2022

We include an excerpt here for convenience:
Digital luminescence images are analyzed using a customized MATLAB (The Mathworks, Inc.) script. The automated workflow is briefly described in the following section. Prior to running the analysis, the user adjusts script parameters to account for the size of the square bounding box region surrounding each sample in the image (estimated number of pixels, width and height). When executed in the MATLAB computational environment (MATLAB version 2021a) the script functions allow the user to select a specific image file for display, and a graphical user interface prompts the user to locate the positions of the upper left sample, the upper right sample, and the lower left sample. The user also enters the total number of sample rows and columns that are in the selected contiguous block of samples. Based on the above information that is supplied by the user, the script functions adjust the image to correct for any rotation of the samples and applies an evenly spaced analysis grid. Each grid region contains luminescence information from a single sample. An automatic threshold is applied to each region to segment the signal region from the local background region within each grid region. The size (in pixels) of the signal region is reported and can be automatically gated to eliminate irregular signals (for example, small false positive speckle noise regions or large false positive regions when no true signal region is present in the grid). Numerical values for mean signal intensity per signal region, mean local background intensity, and mean pixel intensity region minus local background are reported for each analyzed grid region.

## Academic credits

The proof of concept for the use of this Matlab file for High-throughput cellular thermal shift assay (CETSA) measurement is published in a peer-reviewed journal article, which is cited below.

Title:
High-throughput cellular thermal shift assay (CETSA) using acoustic transfer of protein lysates
Authors:
Ashley E. Owens1, Michael J. Iannotti1, Tino W. Sanchez1, Ty Voss1, Abhijeet Kapoor1, Matthew D. Hall1, Juan J. Marugan1, Sam Michael1, Noel Southall1, Mark J. Henderson1* 
Affiliation:
1 National Center for Advancing Translational Sciences, National Institutes of Health, Rockville, Maryland, 20850, USA
*Corresponding Author:  Mark Henderson; email: mark.henderson2@nih.gov

Published in ACS Chemical Biology 2022

Any future projects using this source code or derivative code should cite the above journal article. Full citation information for this journal article will be available at https://www.ncbi.nlm.nih.gov/pubmed
