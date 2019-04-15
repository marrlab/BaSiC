# BaSiC

Matlab code accompanying 

**A BaSiC Tool for Background and Shading Correction of Optical Microscopy Images**

by Tingying Peng, Kurt Thorn, Timm Schroeder, Lichao Wang, Fabian J Theis, Carsten Marr\*, Nassir Navab\*, Nature Communication 8:14836 (2017). [doi: 10.1038/ncomms14836](http://www.nature.com/articles/ncomms14836).

BaSiC is licensed under 

[Creative Commons Attribution-NonCommercial-NoDerivatives 4.0 International Public License](https://creativecommons.org/licenses/by-nc-nd/4.0/legalcode)

It is free for academic use and please contact us for any commercial use.

## Demo

Download demo data examples from [Dropbox](https://www.dropbox.com/s/plznvzdjglrse3h/Demoexamples.zip?dl=0) and run matlab files under example folder.

## ImageJ/Fiji Plugin
BaSiC is also available as a ImageJ/Fiji Plugin.


### Installation instruction

Note: If you do not have Fiji installed on your computer, you can download it from [Fiji website](http://fiji.sc/).


### Install via Fiji Updater

1. Start Fiji and run the updater ("Help->Update Fiji")
2. Select the "Manage Update Sites" button at the bottom-left of the updater window
3. Scroll the list of available update sites to find "BaSiC" (Note: If you cannot find "BaSiC" in the list, select "Add Update Sites", Change the name field from default "New" to "BaSiC", set the URL field to http://sites.imagej.net/BaSiC/)
4. Check the box at the left of "BaSiC"
5. Select "Close" 
6. Select "Apply Changes" 
7. Restart Fiji. BaSiC should appear in the Plugins menu.

From now on, running the Fiji updater will also check for BaSiC updates, and install them if they are available.


### Install manually

Please download [BaSiC Plugin](https://github.com/QSCD/BaSiC/blob/master/BaSiCPlugin.zip) from this repository. 

1. Copy “BaSiC_.jar” to the “$FIJIROOT/plugins” folder of your Fiji/ImageJ installation.
2. Copy all dependent jar files in the "Dependent" folder to your Fiji/ImageJ "$FIJIROOT/jars" directory.


### Troubleshooting

If you get the error message 

"java.lang.NoSuchMethodError: edu.emory.mathcs.utils.ConcurrencyUtils.submit"

make sure that in your Fiji/ImageJ "$FIJIROOT/jars" directory, there is only one version of each jar from the "Dependent" folder. Particularly, delete jtransforms-2.4.jar and replace it with our jtransform.jar.

## Issues
If you have any issues concerning BaSiC, please report them in the [Issues](https://github.com/QSCD/BaSiC/issues) section of this GitHub repository and we will try to find a solution.






