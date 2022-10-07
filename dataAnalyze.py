# "dataAnalyze.py" is a pyROOT script that automates the data analysis for input proton source data saved from
# "FPv56.vi". The previous program "databuild.C" was only meant to generate ROOT TTrees for viewing the data and doing
# basic analysis via TTree->Draw() options.

# Input: ROOT file containing TTree from databuild.C output, time increment of averaging

################### Process #######################

# Step 0: For a more universal analysis, normalize the H2+ count rate time series histogram against the Keller LEO 3 Pressure gauge
#   Keep the H2+ data since we will also simply normalize the proton count rate against the H2+ rate
#   Practically, this step means dividing histograms: ->Divide() function (this is bin by bin division per 2883 of TH1.cxx)


# Step 1: Obtain the average and standard deviations (defined as RMS) of the resulting time series histograms for given
#   time increments. If the total time interval exceeds an integer multiple of the time increment, the last data are excluded.
#   Plot & save the results to a histogram (TH1)

# Step 2: Fit a line of the results to obtain a 'stability = intercept/slope' datum in units of ppm and ppm/hr

####################################################

# Output: Save the resulting analysis (histograms, functions etc.) to a new root file to keep all the data and analysis together.

import sys

# These are imported for quick sum of squares (RMS).
from numpy import dot
from numpy import array

from ROOT import TF1
from ROOT import TH1F
from ROOT import TFile
from ROOT import gStyle

def main(argv):

    ####################### INPUT HANDLING ############################
    
    if len(sys.argv) != 4:
        sys.exit('usage: pyroot dataAnalyze.py ROOTfile INTEGERtimeIncrement timeIncrementUnits')

    if str(argv[0]).find('.root') < 0:
        sys.exit('program only accepts ROOT files')

    print('opening ' + str(argv[0]) + '...')
    
    #   Try to open file given in the argument and exit otherwise
    try:
        infile = TFile.Open(str(argv[0]))
    except:
        sys.exit()
    
    #   These will 'get' the time series histograms from the opened file and return 'nil' otherwise
    Proton      =   infile.Get('ProtonCount')
    H2          =   infile.Get('H2Count')
    KellerBar   =   infile.Get('KellerP_bar')
    
    if not Proton or not H2 or not KellerBar:
        sys.exit('required histograms "ProtonCount", "H2Count", "KellerP_bar" not found')

    if Proton.GetNbinsX() == H2.GetNbinsX():
        NBins = Proton.GetNbinsX() # This gets the total bin numbers
    else:
        sys.exit('expect all time series bins to be the same (check databuild)')

    #   The 'timeIncrement' and 'timeIncrementUnits' (defined in 'usage' check above) are interpreted by selection

    try:
        delta_t = int(argv[1])
    except ValueError:
        sys.exit('second argument (timeIncrement) must be integer.')

    timeUnit = str(argv[2])
    
    ####################### INPUT PROCESSING ############################

    ProtonRatio     =   Ratio(Proton, H2)
    H2Ratio         =   Ratio(H2, KellerBar)
    ProtonRatio2    =   Ratio(Proton, H2Ratio)
    
    averaging_interval = 10
    averaging_units  = 's'
    averaging_unit_convert = interpret_time_unit(averaging_units)
    averaging_time_seconds = averaging_interval * averaging_unit_convert

    timeUnitConvert = interpret_time_unit(timeUnit)
    delta_ts = delta_t * timeUnitConvert

    # Since averaging changes the bin number, this needs to be accounted for when fitting over intervals
    postAvg_delta_ts = delta_ts//averaging_time_seconds

    # The point of averaging is to have control over the bin errors for the fits.
    # The actual bin errors would require propagation from signal through each device.
    averagedHists   = average_intervals(averaging_interval, averaging_units, [ProtonRatio, H2Ratio, ProtonRatio2])
    hists = averagedHists
    
    function = TF1('line1', 'pol1', 1, NBins)
    hists[0].Fit('line1','S','', 0, NBins)
    input('stuff')
    ProtonRatio.Fit('line1','S','', 0, NBins)
    input()
    #hists = [ProtonRatio,H2Ratio,ProtonRatio2]
    fit_results = []
    stabilities = []
   
    for hist in hists:
        fit_results.append(linear_fits(postAvg_delta_ts,'s',hist))
#        print('fitting') 
        stabilities.append(stability(fit_results[-1][0]))
    
    draw_fits(hists[0],fit_results[0][1])
    input("Enter to continue...")

def linear_fits(delta_t, timeUnit, hist):
    # Evidently one would want to use this with a larger time interval than the averaging interval
    
    NBins = hist.GetNbinsX()
    delta_ts = delta_t * interpret_time_unit(timeUnit)

    # The // is floor division (i.e. answer is integer supremum of the result)
    nIncrements = NBins // delta_ts 

    fits = []
    fit_results = []
    BinsExcluded = identify_negative_values([hist])[0] 

    function = TF1('line', 'pol1', 1, NBins)
    


    #def selected_func(): # We would like to reject points from the fit as an alternative to the points exclusion

        # Takes in a TF1 object and excludes bins as requested.
    #    function = TF1('line', 'pol1', 1, NBins)
        #for bin_number in BinsExcluded:
         #   if hist.GetBinContent(bin_number) < 0:
         #       function.RejectPoint(1)

        #return function

    for i in range(nIncrements):
        # Recall bin 0 is underflow. I.e. we need to start at bin 1.
        start_bin   =   (i * delta_ts) + 1
        end_bin     =   start_bin + delta_ts

        # The fit happens here and the resulting objects are appended to a list.
        # 'poln1' is a first order polynomial fit, i.e. a line f(x) = p0 + p1*x
        # The 'S' option specifies to save the fit statistics to the resulting objects,
        # fit_results[k].Chi2() accesses the chi-squared stat, '.Parameter(n)' returns pn, '.ParError(n)' returns the error of pn
        # Our number of degrees of freedom is delta_ts - len([p0,p1]) = (number of bins - number of fit parameters).
        fit_results.append( hist.Fit('line','S Q','', start_bin, end_bin)  )
        # We also want the functions with their domains
        fits.append(hist.GetFunction('line'))

    results = [fit_results,fits]
    return results


def draw_fits(hist,fits):
    # We want some function to plot the histogram with the fits over each time interval @ )--;--'--- 
    gStyle.SetOptStat(0)        
    hist.Draw('HIST')

    nIncrements = len(fits)
    
    for i in range(nIncrements):
        fits[i].Draw('SAME')
    

def stability (fit_results):
    # 'Stability' is defined as intercept/slope such that as slope -> 0 we have arbitrarily high stability.
    # Variability may be defined as slope/intercept such that slope -> 0 means we have zero variability. 
    stabilities         = []
    stabilities_error   = []

    nIncrements = len(fit_results)

    for i in range(nIncrements):
            
        results = fit_results[i]

        slope       =   results.Parameter(1)
        intercept   =   results.Parameter(0)
            
        stability   =   intercept/slope

        # The stabilites in default units is counts/(counts/seconds) = seconds
        stabilities.append(stability)
            
        slope_error     =   results.ParError(1)
        intercept_error =   results.ParError(0)
            
        # Assuming uncorrelated errors, its sufficient to return the second root of the errors added in quadrature
        stability_error = (slope_error**2 + intercept_error**2)**0.5
            
        stabilities_error.append(stability_error)
    
    return [stabilities, stabilities_error]


def Ratio(hist1, hist2):
    # Histograms will get ratioed (Division with renaming)
    histRatio = hist1.Clone()
    
    histRatio.SetName(hist1.GetName() + '_' + hist2.GetName() + 'Ratio')

    histRatio.GetYaxis().SetTitle('Ratio')

    hist1Title = hist1.GetTitle().replace(' Time Evolution', '')
    hist2Title = hist2.GetTitle().replace(' Time Evolution', '')
    
    histRatio.SetTitle('#font[132]{' + hist1Title + '/' + hist2Title  + ' Ratio Time Evolution}')
    
    # Note that x.Divide(TH1F y) returns a Boolean & implements x = x/y, so there is no need for new definitions.
    Bool = histRatio.Divide(hist2)
    
    if not Bool:
        sys.exit('Division ' + hist1.GetName() + '/' + hist2.GetName() + 'failed')

    return histRatio


def identify_negative_values(histlist):
    # Input are assumed to be 1D histograms with the same bin number
    # Want to exclude the '-9' non-readings from the data (-9 being the default failure mode return set by "FPv56.vi")
    # In their place use the mean or mode of the time increment or omit the data from the analysis.
    NHists = len(histlist)
    BinsExcluded = [0.0] * NHists
    for k in range(NHists):

        NBins = histlist[0].GetNbinsX()
        BinsExcluded[k] = []

        if histlist[k].GetNbinsX() != NBins:
            sys.exit('error: differing number of bins')

        for i in range(NBins):
            if histlist[k].GetBinContent(i) < 0:
                BinsExcluded[k].append(i)

    return BinsExcluded


def interpret_time_unit(timeUnit):
    if str(timeUnit) == 's' or str(timeUnit).find('sec') >= 0:
        #print('assuming units of seconds...')
        timeUnitConvert = 1

    elif str(timeUnit) == 'm' or str(timeUnit).find('min') >= 0:
        #print('assuming units of minutes...')
        timeUnitConvert = 60
    else:
        sys.exit('time units must be seconds (s, sec) or minutes (m, min)')
    
    if not timeUnit:
        sys.exit('time units not able to be defined.')

    return timeUnitConvert


def average_intervals(delta_t, timeUnit, histlist):
    # There is a performance cost for looping over one histogram at a time, but there is clarity in calling the function for each histogram
    # individually. The solution is to take in an arbitrary list of histograms and do the averaging for all of them.
    
    timeUnitConvert = interpret_time_unit(timeUnit)
    
    delta_ts = delta_t * timeUnitConvert
    NBins = histlist[0].GetNbinsX() # This gets the total bin numbers

    # NBins default units is seconds, so if the increment delta_t is in minutes we convert it to seconds
    nIncrements = NBins//delta_ts # '//' is floor division (neglect decimal digits)
    MBins = nIncrements * delta_ts # This is our new bin number neglecting bins at one of the ends

    # Simpler to set the bin contents outside of the main loop so we need to store the calculations in lists
    # To do the averaging, calculate the mean and RMS in the loop over bins. To have memory of past values, we use fixed size lists, Values[i].
    Avg = []
    StdDev = []
    Values = [[]]
    AvgHists = []
       
    NHists = len(histlist)
    for i in range(NHists): 
        Avg.append([])
        StdDev.append([])

        if i > 0:
            Values.append([])

        histTitle = histlist[i].GetTitle().replace('Time Evolution','')
        AvgHists.append(
           
            TH1F(histlist[i].GetName() + 'Avg', '#font[132]{' + histTitle + ' Average};#font[132]{Time ' + str(delta_t) + ' ' + str(timeUnit) + '}',
                nIncrements, 0.0, nIncrements)
        )
    # We need to exclude negative values from the averages (these are non-readings '-9's).
    # The function on the right hand side returns a list of excluded bins for each histogram
    # The object on the left hand side is 2D list with max. number of elements [histlist][bins excluded (dynamically)]
    BinExcludeList = identify_negative_values(histlist)

    # Fill the above histograms with the first MBins < NBins of the original histograms
    # The first bin is 1, which is the value shown from [0,1). So the (n+1)th bin is the value in the nth bin.
    for i in range(1, MBins+1):

        for k in range(0, nIncrements):
            
            # This should fill the first 'nIncrements' of these lists with the first nIncrement values of the histograms
            if (i-1) % nIncrements == k:
                for j in range(NHists):

                    # If the bin is in the excluded bin we won't include it...
                    if not (i in BinExcludeList[j]):
                        Values[j].append( histlist[j].GetBinContent(i) )

            # There may be a more universal condition, but the latter statements in each line work given the former inequalities respectively.
            condition1 = (delta_ts < nIncrements and k == delta_ts)
            condition2 = (delta_ts >= nIncrements and Values[0][-1] != 0.0 and Values[0][-1] != 0.0)
            
            if i % delta_ts == 0 and (condition1 or condition2):
                average(Avg,Values)
                std_dev(StdDev,Values)
                
                # Reset the lists after calculating the statistics 
                for j in range(NHists):
                    Values[j].clear()
    
    for k in range(nIncrements):
        for j in range(NHists):
            AvgHists[j].SetBinContent(k+1, Avg[j][k])
            AvgHists[j].SetBinError(k+1, StdDev[j][k]/(delta_ts)**0.5)

    return AvgHists


def average(avglist, valuelist):
    # Take average   
    if len(valuelist) != len(avglist):
        sys.exit('Error, cannot take average.')
    
    for k in range(len(avglist)):
        avglist[k].append(sum(valuelist[k])/len(valuelist[k]))
    #except:
    #    sys.exit('Failed to append average.')


def std_dev(stdevlist, valuelist):
    if len(valuelist) != len(stdevlist):
        sys.exit('Error, cannot take standard deviation.')
    
    nValuelist = [0.0] * len(stdevlist)
    for k in range(len(stdevlist)):
        # Take standard deviation
        # Define vectors to take the dot products
        nValuelist[k] = array(valuelist[k])
                                
        # Numpy arrays will operate on each element with an integer (i.e. [2,3] - 1 = [1,2])
        # Here we want the central second moment (standard deviation) rather than RMS in a strict sense.
        avg = sum(valuelist[k])/len(valuelist[k])
        stdevlist[k].append((dot(nValuelist[k]-avg, nValuelist[k]-avg)/len(valuelist[k]))**0.5)    


if __name__ == "__main__":
    main(sys.argv[1:])
