pro SAUCES_Kepler_min_histplots, directoryname, target, period, periodtext
; 
;  This program will plot each Kepler quarter for the light curves
;  It will then plot a histogram for each quarter for the time
;  differences.
;  There will then be a histogram for the entire Kepler light curve
;
;  Includes output of histograms showing time between minima
;  R. Roettenbacher 7 Aug 2017
;  4 January 2019  Cleaned up to go on github
;
;  This procedure requires being passed the target (name as a string),
;  the period (as a double), and the period (as a string for easy plotting

CLOSE, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10

IF period LT 1. THEN BEGIN
	binsize = 0.025
ENDIF 
IF period GE 1. AND period LT 2.0 THEN BEGIN
	binsize = 0.1
ENDIF
IF period GE 2.0 THEN BEGIN
	binsize = 0.25
ENDIF

permult = 3.

; Read in the light curve
READCOL, directoryname + 'FullLCs_CBV/' + target + '_LC_CBV.txt', mbjd, flux, FORMAT = 'd,d'


; Read in the minimum file
READCOL, directoryname + 'TimeofMinima_CBV/' + target + '_TimeofMinima_CBV.txt', mintime, smoothflux, unsmoothflux, FORMAT = 'd,d,d'


; Read in quarter starts and ends
READCOL, directoryname + 'Kepler_quarter_beginend.txt', quarter, qstart, qend, FORMAT = 'a,d,d'

; Set up the array for the differences between the minima
timediff = dblarr(N_ELEMENTS(mintime) - 1)
	
; Calculating the time difference between minima
FOR j = 0, N_ELEMENTS(mintime) - 2 DO BEGIN
	timediff[j] = mintime[j+1] - mintime[j]
ENDFOR

; Plot the whole light curve and the histogram for its minima
A = FINDGEN(17)*(!PI*2/16.)
USERSYM, COS(A)/5., SIN(A)/5., /FILL

SET_PLOT, 'ps'
DEVICE, FILENAME = directoryname + 'TimeofMinima_CBV_histplots/' + target + '_minhistplots.eps', /LANDSCAPE, /COLOR, DECOMPOSED = 1
!x.style = 1.
;!y.style = 1.
!p.thick = 2.
!x.thick = 2.
!y.thick = 2.
!p.charsize = 1.5
!p.symsize = 2.
!p.font = 2
!p.multi = [0,1,2]


plot, mbjd, flux, psym = 8, title = target + '    ' + periodtext + ' days', ytitle = 'Normalized Flux', xtitle = 'BJD - 2454833', yrange = [min(flux)-0.005, max(flux)+0.005], position = [.1,.75,.9,.9], ystyle = 1

plothist, timediff, bin=binsize, xtitle = 'Time between Minima (days)', ytitle = 'Number of Occurrences',  position = [.1,.1,.9,.6]
	

plot, mbjd, flux, title = target + '    ' + periodtext + ' days', ytitle = 'Normalized Flux', xtitle = 'BJD - 2454833', yrange = [min(flux)-0.005, max(flux)+0.005], position = [.1,.75,.9,.9], ystyle = 1

plothist, timediff, bin=binsize, xtitle = 'Time between Minima (days)', ytitle = 'Number of Occurrences', xrange = [0,permult*period],  position = [.1,.1,.9,.6]

FOR i = 0, N_ELEMENTS(quarter) - 1 DO BEGIN
	test = WHERE(mbjd GT qstart[i] AND mbjd LT qend[i]) 
	IF (N_ELEMENTS(test) GT 1) THEN BEGIN
		indexlc = WHERE(mbjd GT qstart[i] AND mbjd LT qend[i])
		tempmbjd = mbjd[indexlc]
		tempflux = flux[indexlc]

		plot, tempmbjd, tempflux, title = target + '    ' + periodtext + ' days    ' + quarter[i], ytitle = 'Normalized Flux', xtitle = 'BJD - 2454833', yrange = [min(tempflux)-0.005, max(tempflux)+0.005], position = [.1,.75,.9,.9], ystyle = 1

		indexmin = WHERE(mintime GT qstart[i] AND mintime LT qend[i])
		tempmin = mintime[indexmin]
		IF N_ELEMENTS(tempmin) GT 3 THEN BEGIN
			tempdiff = dblarr(N_ELEMENTS(tempmin) - 1)
			; Calculating the time difference between minima
			FOR k = 0, N_ELEMENTS(tempmin) - 2 DO BEGIN
				tempdiff[k] = tempmin[k+1] - tempmin[k]
			ENDFOR
			plothist, tempdiff, bin=binsize, xtitle = 'Time between Minima (days)', ytitle = 'Number of Occurrences', xrange = [0,permult*period],  position = [.1,.1,.9,.6]

		ENDIF ELSE BEGIN
			plothist, timediff, bin=binsize, xtitle = 'Time between Minima (days)', ytitle = 'Number of Occurrences', xrange = [0,permult*period],  position = [.1,.1,.9,.6], /nodata
		ENDELSE
	ENDIF
ENDFOR

CLEANPLOT, /SILENT
DEVICE, /CLOSE
SET_PLOT, 'x'

RETURN
END

