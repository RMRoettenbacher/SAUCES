pro SAUCES_period_find_FFT, directoryname, target, initperiod, rotperiod, signaltype
;
; Uses periodogram.pro to find the most likely period of the light curve
; for the whole LC
;
; R. Roettenbacher 9 October 2017
; 2 June 2018
; 16 November 2018 Cleaned and commented up for github

CLOSE, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10

; Read in the light curve (here, the CBV-corrected light curve)
READCOL, directoryname + '/FullLCs_CBV/' + target + '_LC_CBV.txt', mbjd, flux, FORMAT = 'd,d'

; In Roettenbacher & Vida (2018), we used the McQuillan et al. (2014) periods
; as a starting point for our period search to reduce the computation time.
; Here, we create the range of period values to run our period-search on.
periodstartmin = initperiod - 10.
periodstartmax = initperiod + 10.


; We set a minimum for the rotation period, so each period reported
; will be greater than the time sampling of Kepler long-cadence data.
IF periodstartmin LT 0.125 THEN BEGIN
	periodstartmin = 0.125
ENDIF


; FFT
power = periodogram(mbjd, flux, period_range=[periodstartmin,periodstartmax], npts=10000)

max = WHERE(power[1,*] EQ MAX(power[1,*]))

max = FIX(max)

maxperiod = power[0,max]
maxperiod = DOUBLE(maxperiod[0])

print, 'Rough period:  ', maxperiod

PLOT, power[0,*], power[1,*], xrange = [maxperiod-3., maxperiod+3.]


; We zoom in on the period to get a better estimate.
; This helps to differentiate between periods when there are a cluster of peaks 
; (evidence of differential rotation). 
IF maxperiod GT 3.125 THEN BEGIN
	powerzoom = periodogram(mbjd, flux, period_range=[maxperiod-3.,maxperiod+3.], npts=10000)
ENDIF ELSE BEGIN
	powerzoom = periodogram(mbjd, flux, period_range=[0.125,maxperiod+3.], npts=10000)
ENDELSE

PLOT, powerzoom[0,*], powerzoom[1,*]

maxzoom = WHERE(powerzoom[1,*] EQ MAX(powerzoom[1,*]))
maxzoom = FIX(maxzoom)

maxzoomperiod = powerzoom[0,maxzoom]
maxzoomperiod = DOUBLE(maxzoomperiod[0])

maxzoompower = powerzoom[1,maxzoom]
maxzoompower = DOUBLE(maxzoompower[0])

OPLOT, [maxzoomperiod,maxzoomperiod],[0,maxzoompower*2], linestyle=2

; Find the other strong peaks.
; The strongest peak in a cluster of peaks is not necessarily the 
; one with the shortest period, so we do the following to try to find
; the shortest period in the cluster.  We will assign that as the 
; stellar rotation.  
localmaxes = WHERE(powerzoom[1,*] GT maxzoompower/3.)
localmaxesperiod = powerzoom[0,localmaxes]
localmaxespower = powerzoom[1,localmaxes]


; Account for maxima that are so much higher than the rest of the
; powers.  Note that this indicates a lot of stability in the light curve
; across the observation.  It's likely that this is not the signature 
; of surface evolution (e.g., spots), but maybe pulsations or eclipsing binaries. 
peak = 0.
peakcount = 0.
signaltype = ' '
FOR max = 1, N_ELEMENTS(localmaxes) - 2 DO BEGIN
	IF localmaxespower[max] GT localmaxespower[max-1] AND localmaxespower[max] GT localmaxespower[max+1] THEN BEGIN
		IF peakcount EQ 0 THEN BEGIN
			peak = localmaxesperiod[max] 
			peakcount++
			signaltype = 's'
		ENDIF ELSE BEGIN
			peak = [peak,localmaxesperiod[max]]
			peakcount++
			signaltype = 's'
		ENDELSE
	ENDIF
ENDFOR



IF peakcount EQ 0. THEN BEGIN
	peak = maxzoomperiod
	peakcount++
	signaltype = 'p'
ENDIF


; Assign the rotation period
rotperiod = peak(WHERE(peak EQ MIN(peak)))


PRINT, target, rotperiod


OPENW, 1, directoryname + '/Periods_Types/' + target + '_Period_Type.txt'
PRINTF, 1, target, rotperiod, '  ', signaltype, /IMPLIED_PRINT
CLOSE, 1

OPLOT, [rotperiod,rotperiod],[0,maxzoompower*2]


PLOT, powerzoom[0,*], powerzoom[1,*], xrange = [rotperiod-0.5,rotperiod+0.5]
OPLOT, [maxzoomperiod,maxzoomperiod],[0,maxzoompower*2], linestyle=2
OPLOT, [rotperiod,rotperiod],[0,maxzoompower*2]


firsttitle = rotperiod

IF signaltype EQ 'p' THEN BEGIN
	title =  'Likely not starspots'
ENDIF ELSE BEGIN
	title = 'Likely starspots'
ENDELSE

A = FINDGEN(17)*(!PI*2/16.)
USERSYM, COS(A)/2., SIN(A)/2., /FILL


; Here are some sample plots for checking all works well.
SET_PLOT, 'ps'
DEVICE, FILENAME = directoryname + '/FFT_Plots/' + target + '_fftplot.eps', /LANDSCAPE, /COLOR, DECOMPOSED = 1
!x.style = 1.
!y.style = 1.
!p.thick = 2.
!x.thick = 2.
!y.thick = 2.
!p.charsize = 1.5
!p.symsize = 2.
!p.font = 2
!p.multi = [0,2,2]

PLOT, power[0,*], power[1,*], xrange = [-1.,30], xtitle = 'Period (days)', ytitle = 'Power', title = firsttitle

PLOT, powerzoom[0,*], powerzoom[1,*], xrange = [rotperiod-3,rotperiod+3], xtitle = 'Period (days)', ytitle = 'Power', title = title, /nodata
OPLOT, [maxzoomperiod,maxzoomperiod],[0,maxzoompower*2], linestyle=2, thick = 3, color = 'AAAAAA'x
OPLOT, [rotperiod,rotperiod],[0,maxzoompower*2], thick = 3, color = '777777'x

OPLOT, powerzoom[0,*], powerzoom[1,*]

PLOT, powerzoom[0,*], powerzoom[1,*], xrange = [rotperiod-0.5,rotperiod+0.5], xtitle = 'Period (days)', ytitle = 'Power', title = title, /nodata
OPLOT, [maxzoomperiod,maxzoomperiod],[0,maxzoompower*2], linestyle=2, thick = 3, color = 'AAAAAA'x
OPLOT, [rotperiod,rotperiod],[0,maxzoompower*2], thick = 3, color = '777777'x

OPLOT, powerzoom[0,*], powerzoom[1,*]


index1 = WHERE(ABS(mbjd - 300.) EQ MIN(ABS(mbjd - 300.)))
index2 = WHERE(ABS(mbjd - 335.) EQ MIN(ABS(mbjd - 335.)))
PLOT, mbjd[index1:index2], flux[index1:index2], PSYM = 8, xtitle = 'BJD - 2454833', ytitle = 'Normalized Flux'




CLEANPLOT, /SILENT
DEVICE, /CLOSE
SET_PLOT, 'x'


RETURN
END