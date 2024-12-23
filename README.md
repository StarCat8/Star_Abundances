Sbundance is a jupyter notebook meant to perform automatic spectral analysis of K metal poor stars by using the equivalent width method. 
 Its main features are the stellar parameter space exploration in order to obtain the right stellar atmosphere to apply to the analysis. In order to do so it integrates three other codes, namely abundance, a SPECTRUM routine by R. Gray, ARES by S. Sousa and PyKMOD, by the github user kolecki. Then it uses a gradient-descent (with a random search spin!) like approach to find the absolute minima of the parameter space.

 The way it's meant to be used is the following:
 ---
 Download all the files needed:
 - the latest version of sbundance
 - the pykmod by kolecki
 - the spectrum277c folder*
 - ARES
*as for the spectrum folder while all the main features of spectrum programs the codes are untouched a little bit of them was rewritten mostly, if not solely, to manage the input and the outputs in an automatic fashion. Aside from abundance.c and spectrum.c (which now have their own parameter file to handle inputs) also spaux.c had minor modification, namely the ggets() function. So while you could download the original spectrum files and them make your own minor changes, it is way easier for you to just grab the properly modified files from here.

If everything is set up and installed you open your jupyter notebook and open sbundance.
---
- Before running any cells you should set your computer's paths in the first cell. From path_to_pythonEnv=... to path_to_spectrum=... you should enter the paths present on YOUR machine. The names should be self explanatory, there are also comments in the notebook to help ;). 
- Run the first cell which contains all the basic functions needed.
- Run the third block (for the second block look below at PERFORMING A CONTINUOUS SEARCH***)
  

 The third block defines a Stella class. 
 ---
 This eats a fits file (equipped with a spectrum) and create a Stella object containing a variety of useful info and functions. When a star=Stella('file.fit') is defined a folder in the output directory is created with the name of the star contained in the fits file. You have the basic atmospherical parameters as star.t_eff, star.logg, star.metal, star.v_m.

 As of today (23/12/24) the code is meant to use HARPS-N fits data, this means that you may need to adjust some functions in order to account for different headers name in different fits file (e.g. for the HARPS-N fits the name of the star is under the header "HIERARCH TNG OBS TARG NAME" and is not always like that).
 HARPS-N automatically provide a radial velocity estimate, which is automatically loaded in the star.v_r position, which you can set manually anyway. 

Available methods are:
--
.spectrum() - this display the spectral data contained in the fits file

.dopcor() - this corrects for doppler shift, using the value of star.v_r

.creaModel() - creates an interpolated atmosphere via pykmod and sets up param file for abundance

.synth_creaModel() - creates an interpolated atmosphere via pykmod and sets up param file for spectrum, create the stdatom.dat file, which contains the chemical composition of the star.

.searchParam2() - this method start the atmospherical parameter search. It decomposes the search in 4 searches managed contemporarely with parallel subprocessing, hence it requires at least 4 cores to be ran. .searchParam() and .searchParam1() follows similar philosphy but for 1 and 2 cores respectively, but as of today (23/12/24) the .searchParam2() is the most refined, and 4 cores cpus are pretty much ubiquitous.

.normalize_spectrum() - this method normalizes the spectrum. In order to do so the following algorithm was applied. Divide the full spectrum in chunks of N_avepoint, say 500 points. The continuum should in principle be around the maximum value of the flux, hence discarding the absorption lines. This was done by iteratively setting all values of intensity below a threshold (given by the intensity average in the chunk) to the value of the threshold itself.  In order to account for upward spikes (eg cosmic rays) values significantly above the average are also discarded. This leaves with value fluctuating around the continuum. Of these points the max are chosen, and then averaged with a moving averaging window of customizable dimension (default is 15). Eventually the point are slightly shifted down by a value comparable to half of the variance. This fairly normalizes the spectrum, namely peaks and troughs are excluded by the averaging, but surely more work is needed

.synth() - this method create a syntetic spectrum. It accepts a plethora of arguments, from initial to final wavelenght, to the integration step, which cpu to use, which value of vsini to use to rotationally broaden the lines. All the elements from H to U can be changed (make sure to set their abundances according to the prescription reported in the spectrum documentation, i.e. the metal abundances scales with the [M/H] value). The defaul values are solar abundances. To set an abundance just pass the element as an argument with its abundance, eg star.synth(..., Na=-5.32, ...). This method produces a synthetic spectrum saved in the star.folder directory (ie the "path_to_output + star.nome" folder) with .synstar extension. If "broad" is also a suffix the spectrum was broadened with avsini.

.avsini() - this method just calls the avsini code in spectrum to broaden the lines to account for rotational effects.

.stampaStella() - this method prints neatly the results from the searchParam(), printing the atmospheric parameters and the uncertainties of the results in a text file.


*** More on .searchParam2() ***
-------
 To give a brief sense on how the search is carried out. 
The search is done via the equivalent width method. This is considered enforced if some functions are 0. The function employed are 4, and all of them are dependant on X = (T, logg, [M/H], vturb). They are s(X), F(X), Y(X), S(X). The s(X) gives the slope of the potential energy balance, F(X) is the discrepancy of iron abundance as derived from FeI and FeII lines, Y(X) is a consistency check as the [M/H] is assumed to be given by [Fe/H], so it is the difference between these two values, S(X) is the slope of iron abundance against EWR, where EWR is the reduced equivalent width. If these 4 functions are 0 all at once the we found the correct X! This is no easy feat as it's impossible to set these exactly to 0. In order to find the best values I devised a metric, M, which maps these four functions in a single real value to discern between atmospheres. This value is 0 iff all the functions are 0. The atmosphere to evaluate are chosen as follows.
 We start with 4 random seed distributed in four different region of the parameter space (and chosen to be as physical as possible, eg to high temp corresponds higher logg), these creates 4 different atmospheres (step 0 atmospheres). Then sbundance create the step 1 atmospheres, which updates the atmospherical parameters (X) according to a previously calibrated matrix. From the 2nd step on the atmospherical parameters are updated according to a matrix created during runtime which is basically the Jacobian the funcions s, F, Y, S expressed as a Ist order taylor expansion in the X variable. Inverting the matrix gives the new X
that would set s, F, Y, S to 0 if the I order was exact (but it's not and it is an iterative process anyway). If more than 4 cores are given (I usually run this on a 12 cores cpu) 2 more atmospheres per seed are computed, and these are randomly distributed in a 4d sphere around the new X updated via the matrix, with a radius varying according to the distance from the previous atmosphere (ie if the guessed temperature is changed by 50K from the previous one, a radius of 25K is applied). Iterating this process has shown to consistently reduce the value of M and that usually two, and at least one, search branches end up in the global minimum. The metric M is devised as such that when it gets around the value of 1 we can be pretty sure to be in the neighborhood of the global minimum, and when below 0.1-0.05 X is as close to X_true as one can reasonably hope (ie, the change in effective temperature would be less than 1K, in logg less than 0.01 dex as well for [M/H] and in vturb < 0.01). Actually the vturb direction in parameter space could be probably dropped off and just use an analytical formula for it as the Mashonkina one.
 This method also allows the user to insert the contribution of NLTE effects. This comes with some caveats. The biggest of them is how trustworthy the correction I used are (Amarsi 2016 corrections). These one now implemented are computed on a grids with pretty wide meshes, and also the parameter space sbundance is meant to probe is not entirely superimposed to the grid. I extended the grid to continuum values via a first order interpolation, but I'm not entirely sure of the results. Anyway one can easily turn on or off the NLTE corrections via the flag NLTE (=0 for LTE and =1 for NLTE). 
 Incoming changes in .searchParam2():

 
It is probably more sensible to drop the random search part until we can be sure we're near the global minimum. This way one can allow for more seeds, say, 12 seeds to run separately and only when one gets comfortably close to the global minimum (ie M<1/0.5) redirect all the cores to a random-ish search in the global minimum neighborhood. This would help in avoiding getting stuck in local minima.

On a cold rainy night..
----
 Say a friend of friends, Peppe Panella, is a sleepwalker and lives near the woods which is notoriously a place where "hic sunt leones". It's a cold rainy sunday night and you need to save him to make your friend, Ugo Bugo, happy. Living in a capitalistic world the chances of survival of Peppe crucially depend on the amount of money YOU have, thus how many private investigators you can afford to hire (the cores). The plan is simple, we know where he lives, he cannot be that far away (we know it's a red giant). Unfortunately the rain washed away Peppe's tracks in the open fields in front of the woods so everyone enters the woods from different places (someone even dropped off an helicopter in the middle of the woods!). Eventually detectives start to find tracks below the thick leaves and start following them (gradient descent), but the rain washed away some of those leaving couple of rescuers running in circle for 10 minutes. The detective who dropped off the helicopter though seems on a good run and following tracks he eventually find a pair of socks! "Surely! They're his! He must be close as one does not simply walk without a sock!" - he shouts, everyone came running and start looking for Peppe. Eventually they find him. Except it's not actually him, not necessarily. You see one important point I is that you don't really know Peppe Panella. Moreover there's like a multitude of Peppe Panella's look alike, all living near the woods and all being sleepwalker. So what you actually do is finding the Peppe Panella that most resemble a polaroid picture of Peppe. This is done thanks to a device called the Peppèmetro. Well this Peppe you found has a value of 0.1 in the Peppèmetro, you notice that the heights corresponds and his haircut also, though a mole on his cheeks betray some differences...
 
 it's late in the night now, I guess Ugo Bugo won't even notice...

 
///OLD VERSIONS
 The analysis comes in two flavour, a grid search on a discretized parameter space (making use of pre-computed ATLAS9 stellar atmospheres provided by F. Castelli) and a search in the continuous parameter space, created by interpolation by the means of the code PyKMOD.
  To code is meant to work with .fits 1-d spectra where the doppler shift is already removed (future works will provide an automatic way to calculate and remove doppler shift).
It is structured in a versatile folder structure, hence the user is required to set the folders' path inside the first block of the notebook (e.g. path_to_output), comments in the code are meant to explain which folder is which.
 The core of the code is the first block, where all the functions are defined in such a way that minimal input are required to the user to obtain the final results.
------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
PERFORMING A GRID SEARCH (not available anymore)

 To perform a full grid search the first function to call is searchPatch() which requires you to only specify the .fits file name. Then, to find out which atmosphere is more likely to be the right one, the function to call is tra_le_atmosfere() or folderLog(). This creates a log file containing all the useful results to evaluate to goodness of the applied atmosphere.
  Lastly call the final_atmo() function to automatically create a "final_result.txt" file containing the parameters of the best atmosphere that was found.
Other methods of grid search are included, but strongly NOT recommended. The cell ###ITERATIVE SEARCH### does an iterative search, evaluating the results every time a new atmosphere is guessed. This has an high chance of getting stuck in local minima.

PERFORMING A CONTINUOUS SEARCH***

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
OLD VERSION///


REQUIREMENTS
 Run pip install -r requirements.txt. (Make sure to have the correct requirement.txt file)

The code makes use of the following libraries:

matplotlib
numpy
astropy
symps
scipy == 13.1 #pykmod should be changed to make it compatible with newer releases
os
re
time
urllib3

The urllib3 library is used only in the scraping functions, you can comment it if you're not going to use it.
