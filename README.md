Sbundance is a jupyter notebook meant to perform automatic spectral analysis of K metal poor stars by using the equivalent width method. 
 It's main features are the stellar parameter space exploration in order to obtain the right stellar atmosphere to apply to the analysis. In order to do so it integrates three other codes, namely abundance, a SPECTRUM routine by R. Gray, ARES by S. Sousa and PyKMOD, by the github user kolecki. Then it uses a gradient-descent like approach to find the absolute minima of the parameter space.

 The analysis comes in two flavour, a grid search on a discretized parameter space (making use of pre-computed ATLAS9 stellar atmospheres provided by F. Castelli) and a search in the continuous parameter space, created by interpolation by the means of the code PyKMOD.
  To code is meant to work with .fits 1-d spectra where the doppler shift is already removed (future works will provide an automatic way to calculate and remove doppler shift).
It is structured in a versatile folder structure, hence the user is required to set the folders' path inside the first block of the notebook (e.g. path_to_output), comments in the code are meant to explain which folder is which.
 The core of the code is the first block, where all the functions are defined in such a way that minimal input are required to the user to obtain the final results.
------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
PERFORMING A GRID SEARCH

 To perform a full grid search the first function to call is searchPatch() which requires you to only specify the .fits file name. Then, to find out which atmosphere is more likely to be the right one, the function to call is tra_le_atmosfere() or folderLog(). This creates a log file containing all the useful results to evaluate to goodness of the applied atmosphere.
  Lastly call the final_atmo() function to automatically create a "final_result.txt" file containing the parameters of the best atmosphere that was found.
Other methods of grid search are included, but strongly NOT recommended. The cell ###ITERATIVE SEARCH### does an iterative search, evaluating the results every time a new atmosphere is guessed. This has an high chance of getting stuck in local minima.

PERFORMING A CONTINUOUS SEARCH

To explore the continuous parameter space it is sufficient to launch the *continuous search block*, making sure to provide the correct name of the star to be analysed and specifing an output folder (to not be confuse with main output folder pointed at by the path_to_output variable).
This makes use of the function avvicina_interpol3(). This function is built in such a way that NLTE corrections can be included in the exploration during run-time, though a calibration is required. Future works will include line-by-line NLTE corrections, to do so a simple scrape function is built in createCorrectGrid() which makes use of the INSPECT database. (careful, since a lot of call will be made to the site in all probability the connection will be refused after quite some time, basically being a DDoS attack. Try asking maybe ;)   )
 Similarly to the grid search a log file containing the results for all the analysed atmospheres is created and then the best atmospheres are collected in the "final_result.txt" file.

To code is meant to be versatile and basically every variable can be manually adjusted to accomodate the users need. It is possible to set a range of accepted EW by changing the values by passing the arguments eqw_low and eqw_high to the preparaAbundance() function, 10 and 150 respectively are default.
 It has built in a way to launch multiple sub-processes of the abundance code, this was made to pave the way to future parallelizations. The code has an original atmosphere interpolator, just call the function creaInterpolazione() specifing the temperature, log(g) and metallicity, although the PyKMOD interpolator works much better (some refinement to the function are still required).



We suggest for both the grid and the continuous search to provide a line-list containing only FeI and FeII lines, these lines being the only ones used to determine the stellar parameter, since the analysis of many lines may results in long computation time by abundance. 100 lines are usually computed in around 20 to 40 seconds on my machine (which is sadly a very poor one).
-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
DERIVING STELLAR ABUNDANCES

The complete iter to obtain stellar abundances is the following:
First you have to perform either a grid search or a continuous search to find the correct atmospheric parameters using the only FeI and FeII line-list. Then, once the atmosphere is obtained, launch the crea_ATMinterpol() using a line-list with all the required elements by passing the argument lineList to the preparaAbundance() function. In the folder uploaded here on github you'll find the line-list. The one for FeI and FeII is called Iron_lines, while another line-list is molte_linee containing alpha and iron peak elements.
 The output file is created. To read it you can use the function simpleAnalysisSingle() which requires you to provide the path and the name of the output (though default values are set), moreover you can specify which elements to check out by specifing the species code through the argument elemento, 26.0 (FeI) is default.
 To make sure the results are somewhat trustworty make sure to check out the values for "s(FeI)" and "FeI-FeII" in the log files or equivalently the last two columns of the final_result.txt files are sufficiently small, <=0.002 for the "FeI-FeII" and <0.001 for "s(FeI)". In theory the smallest the better.

The crea_ATMinterpol() function uses interpolated atmospheres, this way all values of atmospheric parameters are allowed. One may want to use the gridded atmospheres (they are actually more trustworthy, though the mesh of the grid may be too wide for use being separated by 250 kelvin in temperature and .5 dex in both log(g) and metallicity). If you want to just run the analysis using a gridded atmosphere you have to copy the atmosphere file inside the "spectrum277c" folder with the name "ModelE". Then run the SlanciaAbundance() function. This function takes 1 argument, namely the number of sub-processes running abundance. 1 is default. To make use of multiple sub-processes you have to copy the atmospheres of interest in the spectrum277c_0, spectrum277c_1 and spectrum277c_2 folders if 2, 3 or 4 sub-processes are required respectively. (e.g. one atmosphere in spectrum277c and one in spectrum277c_0 if two sub-processes are required, one in each four folder if 4 sub-processes). 
-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------


REQUIREMENTS

The code makes use of the following libraries:

matplotlib
numpy
astropy
scipy
os
re
time
urllib3

The urllib3 library is used only in the scraping functions, you can comment it if you're not going to use it.
