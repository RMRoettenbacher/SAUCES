; pro SAUCES_Kepler_flares_compare_test.pro
; 
; Reads in the timing of the flares detected by FLATW'RM and compares them to the minima 
; detected by the SAUCES pipeline.
;
; R. Roettenbacher 7 May 2018
; 4 January 2019 Cleaned up for github

directoryname = '/Your/Directory/Name/'

FILE_MKDIR, directoryname + 'Flare_Compare'
FILE_MKDIR, directoryname + 'Flare_Compare_Hist'
FILE_MKDIR, directoryname + 'All_Flare_Compare'
FILE_MKDIR, directoryname + 'All_Flare_Compare_Hist'

A = FINDGEN(17)*(!PI*2/16.)
USERSYM, COS(A)/1., SIN(A)/1., /FILL, COLOR = '000000'x

; Your target list here!
READCOL, directoryname + 'Kepler_flare_target_list.txt', kepID, FORMAT = 'a', /SILENT

; Your period list here!  
READCOL, directoryname + 'Kepler_flare_period_list.txt', FFTper, FORMAT = 'd', /SILENT

; Your mass list here!  
READCOL, directoryname + 'Kepler_flare_mass_list.txt', kepmass, FORMAT = 'd', /SILENT

; Your list of star types (output from SAUCES script)
READCOL, directoryname + 'Kepler_flare_FFT_kepid_period_typeofstar_list.txt', name, startype, FORMAT = 'a,a', /SILENT

; Use these for different cuts in flare strength 
flarelowerlimit = [0.01,0.05]
flareupperlimit = [10.0,10.0]
flarelist = ['001-10','005-10']

; If a flare happens this close to a minimum, tells code to look at adjacent rotations 
jdfix = [0.40]
jdlist = ['040']

endjdfix = [0.40]
endjdlist = ['040']

bins = 100
					
bottomleftx = (FINDGEN(bins)-50.)/100. 
bottomrightx = (FINDGEN(bins)-50.+1.)/100.
topleftx = bottomleftx
toprightx = bottomrightx
	
bottomlefty = DBLARR(bins)
bottomrighty = DBLARR(bins)
toplefty = DBLARR(bins)
toprighty = DBLARR(bins)

numflares = DBLARR(3,bins)

spots = 0
flarecount = 0
actualspotcount = 0

previntegral = 0.
FOR b = 0, N_ELEMENTS(jdfix) - 1 DO BEGIN
FOR c = 0, N_ELEMENTS(endjdfix) - 1 DO BEGIN



	PRINT, jdfix[b], endjdfix[c]
	allflarefluxes = 0.
	allphasediffs = 0.
	keepflarefluxes = 0.
	keepphasediffs = 0.
	keepflarefluxes5 = 0.
	keepphasediffs5 = 0.

	FOR i = 0, N_ELEMENTS(kepID) -1 DO BEGIN
		READCOL, directoryname + 'Periods_Types/' + kepID[i] + '_Period_Type.txt', values, FORMAT = 'a', /SILENT	
		target = values[0]
		rotperiod = values[1]
		FFTper[i] = rotperiod
		signaltype = values[3]

				
		IF signaltype EQ 's' AND startype[i] EQ 's' THEN BEGIN ;  AND FFTper[i] GT 5.0  AND kepmass[i] LE 0.5 THEN BEGIN ; 

			spots++
			READCOL, directoryname + 'TimeofMinima_CBV/' + kepID[i] + '_TimeofMinima_CBV.txt', minjd, minfluxsmooth, minfluxunsmooth, FORMAT ='d,d,d', /SILENT

			; You can change this depending on how many data points you required FLATWRM to use to identify flares
			READCOL, directoryname + 'FLATWRM_flarepts3_output/' + kepID[i] + '_flarepts3_output.txt', flarestart, flareend, flaremax, flareflux, flareint, fitamp, fitfwhm, fitstart, fitend, fitmax, fitint, FORMAT = 'd,d,d,d,d,d,d,d,d,d,d', /SILENT

			IF previntegral NE 0. THEN BEGIN
				IF flareint[0] EQ previntegral THEN BEGIN
					FOR f = 0, N_ELEMENTS(flareflux) -1 DO BEGIN
						flareflux[f] = 0.
					ENDFOR
				ENDIF
			ENDIF

			previntegral = flareint[0]
			
			READCOL, directoryname + 'FullLCs_CBV/' + kepID[i] + '_LC_CBV.txt', lcdate, lcflux, FORMAT = 'd,d', /SILENT

			READCOL, directoryname + 'Missingdata/' + kepID[i] + '_MissingData.txt', missingstart, missingend, FORMAT = 'd,d', /SILENT

			IF FFTper[i] GT 1.0 THEN BEGIN
				phasediff = DBLARR(N_ELEMENTS(flaremax))
				jddiff = DBLARR(N_ELEMENTS(flaremax))
				endphasediff = DBLARR(N_ELEMENTS(flaremax))
				endjddiff = DBLARR(N_ELEMENTS(flaremax))
				closestmin = DBLARR(N_ELEMENTS(flaremax))
		
				FOR j = 0, N_ELEMENTS(flaremax) - 1 DO BEGIN
					
						index = WHERE(ABS(minjd-flaremax[j]) EQ min(ABS(minjd-flaremax[j])) )
				
						closestmin[j] = minjd[index[0]]

						jddiff[j] = (closestmin[j] - flaremax[j])
						phasediff[j] = jddiff[j]/FFTper[i]
	
	
						endindex = WHERE(ABS(minjd-flareend[j]) EQ min(ABS(minjd-flareend[j])) )
	
						endjddiff[j] = (minjd[endindex[0]] - flareend[j])
						endphasediff[j] = endjddiff[j]/FFTper[i]
	
				
						IF phasediff[j] GT 0.5 AND phasediff[j] LE 1.0 THEN BEGIN
							phasediff[j] = 1. - phasediff[j]
							; If there is a phasediff greater than 1, then that will because the there
							; is a significant gap in the data.  We will ignore any flare that 
							; is further away from a minimum than one rotation.
						ENDIF
						IF phasediff[j] LT -0.5 AND phasediff[j] GE -1.0 THEN BEGIN
							phasediff[j] = 1. + phasediff[j]
						ENDIF
						IF endphasediff[j] GT 0.5 AND endphasediff[j] LE 1.0 THEN BEGIN
							endphasediff[j] = 1. - endphasediff[j]
						ENDIF	
	


						beforeindex = WHERE(ABS(minjd-flaremax[j] + FFTper[i]) EQ min(ABS(minjd-flaremax[j] + FFTper[i])) )
						afterindex = WHERE(ABS(minjd-flaremax[j] - FFTper[i]) EQ min(ABS(minjd-flaremax[j] - FFTper[i])) )
						beforeclosestmin = minjd[beforeindex[0]]
						afterclosestmin = minjd[afterindex[0]]

						spotindex = WHERE(minjd GT beforeclosestmin AND minjd LT afterclosestmin) ;note that the start and end minima aren't included
						; spotindex should EQ 1 if there is one spot
						; spotindex should be GT 1 if there are more spots		
				
						IF N_ELEMENTS(spotindex) NE 1 THEN BEGIN
							phasediff[j] = 100.
							endphasediff[j] = 100.
						ENDIF


						; Get rid of the flares that happen too early or too late in the light curve 
						; (less than half a period after the start or before the end), as these
						; may give erroneous values.

						IF (flaremax[j] LT (lcdate[0] + FFTper[i]/2.)) OR (flaremax[j] GT (lcdate[N_ELEMENTS(lcdate) - 1] - FFTper[i]/2.)) THEN BEGIN
							phasediff[j] = 100.
							endphasediff[j] = 100.
						ENDIF 

						IF (flareend[j] LT (lcdate[0] + FFTper[i]/2.)) OR (flareend[j] GT (lcdate[N_ELEMENTS(lcdate) - 1] - FFTper[i]/2.)) THEN BEGIN
							phasediff[j] = 100.
							endphasediff[j] = 100.
						ENDIF 
					
						; Get rid of the events that happen within two data points of a gap
						; There's a non-negligible chance that FLATW'RM marked these as flares.
						missingcheck = WHERE(flaremax[j] LT missingend + 0.041 OR flareend[j] GT missingstart - 0.021)
						IF missingcheck[0] GT 0. THEN BEGIN
							phasediff[j] = 100.
							endphasediff[j] = 100.
						ENDIF


						missingmincheck = WHERE(closestmin[j] LT missingend + 0.1 AND closestmin[j] GT missingstart - 0.1)
						closetogap = WHERE(flarestart[j] GT missingend AND flarestart[j] LT missingend + 0.1)
						; starts close to, but not in gap...this will send the code to calculating the minima in the next and
						; previous rotations when checked below
	

						
		
						; Some of the flares fall in the dips and cause issues with where the
						; minimum is supposed to be.  Account for that here.  Also, some of the larger 
						; flares distort the shapes of the spot minima significantly.  When looking
						; at either the max/start time or the end time, consider a phase difference
						; range that could be affected and apply the below fix
			
				
						flarecheck1 = WHERE(flaremax GT (closestmin[j] - jdfix[b]) AND (flaremax LT (closestmin[j] + endjdfix[c]))) 
						flarecheck2 = WHERE(flareend GT (closestmin[j] - jdfix[b]) AND (flareend LT (closestmin[j] + endjdfix[c])))
						
						IF flarecheck1[0] GE 0. OR flarecheck2[0] GE 0. OR missingmincheck[0] GE 0. OR closetogap[0] GE 0. THEN BEGIN
							; Check where the minima are before and after the flare are in 
							; order to decide if the flare is skewing the location of the spot
							beforeindex = WHERE(ABS(minjd-flaremax[j] + FFTper[i]) EQ min(ABS(minjd-flaremax[j] + FFTper[i])) )
							afterindex = WHERE(ABS(minjd-flaremax[j] - FFTper[i]) EQ min(ABS(minjd-flaremax[j] - FFTper[i])) )


							; If there is more than one spot on the star, the previous or next 
							; minima won't be one rotation away.  This accounts for that.
			
							beforeclosestmin = minjd[beforeindex[0]]
							afterclosestmin = minjd[afterindex[0]]



	
							IF ((afterclosestmin - beforeclosestmin)/2. GT 0.9*FFTper[i]) AND ((afterclosestmin - beforeclosestmin)/2. LT 1.1*FFTper[i]) THEN BEGIN
								
								; This means that the flare is near the center of a minimum, and
								; there's nothing wrong with the location of the next and previous 	
								; minima (none are missing).  This will then use minima around the flare
											
								beforejddiff = ((beforeclosestmin + FFTper[i]) - flaremax[j])
								afterjddiff = ((afterclosestmin - FFTper[i]) - flaremax[j])

								beforephasediff = beforejddiff/FFTper[i]
								afterphasediff = afterjddiff/FFTper[i]



								closestmin[j] = MEAN([beforeclosestmin,afterclosestmin])
								jddiff[j] = MEAN([beforejddiff,afterjddiff])
								phasediff[j] = MEAN([beforephasediff,afterphasediff])


						IF phasediff[j] GT 0.5 AND phasediff[j] LE 1.0 THEN BEGIN
							phasediff[j] = 1. - phasediff[j]
							; If there is a phasediff greater than 1, then that will because the there
							; is a significant gap in the data.  We will ignore any flare that 
							; is further away from a minimum than one rotation.
						ENDIF
						IF phasediff[j] LT -0.5 AND phasediff[j] GE -1.0 THEN BEGIN
							phasediff[j] = 1. + phasediff[j]
						ENDIF


		
								updateminindex = FLOOR(MEAN([beforeindex,afterindex]))
								minjd[updateminindex] = closestmin[j]


								; If there's a flare in the next or last minima, move one more over...
					
								; First, check that this isn't the first or last flare...
								IF (j GE 1) AND (j LT N_ELEMENTS(flaremax) - 2) THEN BEGIN
							
								; Check to see if there is a flare near either the minimum before or the minimum after
				
									beforeflarecheck1 = WHERE(flaremax GT (beforeclosestmin - jdfix[b]) AND (flaremax LT (beforeclosestmin + endjdfix[c])))
									beforeflarecheck2 = WHERE(flareend GT (beforeclosestmin - jdfix[b]) AND (flareend LT (beforeclosestmin + endjdfix[c])))
									afterflarecheck1 = WHERE(flaremax GT (afterclosestmin - jdfix[b]) AND (flaremax LT (afterclosestmin + endjdfix[c])))
									afterflarecheck2 = WHERE(flareend GT (afterclosestmin - jdfix[b]) AND (flareend LT (afterclosestmin + endjdfix[c])))
									missingbeforemincheck = WHERE(beforeclosestmin LT missingend + 0.1 AND beforeclosestmin GT missingstart - 0.1) 
									missingaftermincheck = WHERE(afterclosestmin LT missingend + 0.1 AND afterclosestmin GT missingstart - 0.1)

									IF beforeflarecheck1[0] GE 0. OR beforeflarecheck2[0] GE 0. OR afterflarecheck1[0] GE 0. OR afterflarecheck2[0] GE 0. OR missingbeforemincheck[0] GE 0. OR missingaftermincheck[0] GE 0. THEN BEGIN
				 
										beforeindex = WHERE(ABS(minjd-flaremax[j] + 2.*FFTper[i]) EQ min(ABS(minjd-flaremax[j] + 2.*FFTper[i])) )
										afterindex = WHERE(ABS(minjd-flaremax[j] - 2.*FFTper[i]) EQ min(ABS(minjd-flaremax[j] - 2.*FFTper[i])) )
		
										; If there is more than one spot on the star, the previous or next 
										; minima won't be one rotation away.  This accounts for that.
		
										beforeclosestmin = minjd[beforeindex[0]]
										afterclosestmin = minjd[afterindex[0]]
				
	
										IF ((afterclosestmin - beforeclosestmin)/4. GT 0.9*FFTper[i]) AND ((afterclosestmin - beforeclosestmin)/4. LT 1.1*FFTper[i]) THEN BEGIN
										
		
											beforejddiff = ((beforeclosestmin + 2.*FFTper[i]) - flaremax[j])
											afterjddiff = ((afterclosestmin - 2.*FFTper[i]) - flaremax[j])
														

											beforephasediff = beforejddiff/FFTper[i]
											afterphasediff = afterjddiff/FFTper[i]
		
		
											closestmin[j] = MEAN([beforeclosestmin,afterclosestmin])
											jddiff[j] = MEAN([beforejddiff,afterjddiff])
											phasediff[j] = MEAN([beforephasediff,afterphasediff])
		


						IF phasediff[j] GT 0.5 AND phasediff[j] LE 1.0 THEN BEGIN
							phasediff[j] = 1. - phasediff[j]
							; If there is a phasediff greater than 1, then that will because the there
							; is a significant gap in the data.  We will ignore any flare that 
							; is further away from a minimum than one rotation.
						ENDIF
						IF phasediff[j] LT -0.5 AND phasediff[j] GE -1.0 THEN BEGIN
							phasediff[j] = 1. + phasediff[j]
						ENDIF


											updateminindex = FLOOR(MEAN([beforeindex,afterindex]))
											minjd[updateminindex] = closestmin[j]

										ENDIF
	
									ENDIF
								
								ENDIF
							ENDIF 
						ENDIF	
							
				ENDFOR
			ENDIF
			IF N_ELEMENTS(flaremax) GT 1 THEN BEGIN
				IF allflarefluxes[0] EQ 0. THEN BEGIN
					tempindex = WHERE(phasediff LE 0.5 AND flareflux GE 0.01)

					IF tempindex[0] GE 0 THEN BEGIN
						actualspotcount++
						allflarefluxes = flareflux(WHERE(phasediff LE 0.5))
						allphasediffs = phasediff(WHERE(phasediff LE 0.5))
						keepflarefluxes = flareflux(WHERE(phasediff LE 0.5 AND flareflux GE 0.01))
						keepphasediffs = phasediff(WHERE(phasediff LE 0.5 AND flareflux GE 0.01))
						keepflarefluxes5 = flareflux(WHERE(phasediff LE 0.5 AND flareflux GE 0.05))
						keepphasediffs5 = phasediff(WHERE(phasediff LE 0.5 AND flareflux GE 0.05))
						PRINT, kepid[i], N_ELEMENTS(WHERE(phasediff LE 0.5 AND flareflux GE 0.01))
	
					ENDIF
				ENDIF ELSE BEGIN
					tempindex = WHERE(phasediff LE 0.5 AND flareflux GE 0.01)

					IF tempindex[0] GE 0 THEN BEGIN
						actualspotcount++
						allflarefluxes = [allflarefluxes, flareflux(WHERE(phasediff LE 0.5))]
						allphasediffs = [allphasediffs, phasediff(WHERE(phasediff LE 0.5))]
						keepflarefluxes = [keepflarefluxes, flareflux(WHERE(phasediff LE 0.5 AND flareflux GE 0.01))]
						keepphasediffs = [keepphasediffs, phasediff(WHERE(phasediff LE 0.5 AND flareflux GE 0.01))]
						keepflarefluxes5 = [keepflarefluxes5, flareflux(WHERE(phasediff LE 0.5 AND flareflux GE 0.05))]
						keepphasediffs5 = [keepphasediffs5, phasediff(WHERE(phasediff LE 0.5 AND flareflux GE 0.05))]
						PRINT, kepid[i], N_ELEMENTS(WHERE(phasediff LE 0.5 AND flareflux GE 0.01))
					ENDIF
				ENDELSE 

				FOR a = 0, N_ELEMENTS(flarelowerlimit) - 1 DO BEGIN
					IF flaremax[0] NE 0. AND N_ELEMENTS(WHERE(flareflux GE flarelowerlimit[a] AND flareflux LT flareupperlimit[a] AND phasediff LT 10.)) GT 1 THEN BEGIN
		
						SET_PLOT, 'ps'
						DEVICE, FILENAME = directoryname + 'Flare_Compare/' + kepID[i] + '_flare_flux' + flarelist[a] + '_jdfix' + jdlist[b] + '_endjdfix' + endjdlist[c] + '_2minmove_100bins_bw_plot.eps', /LANDSCAPE, /COLOR, DECOMPOSED = 1
						!x.style = 1.
						!y.style = 1.
						!p.thick = 2.
						!x.thick = 2.
						!y.thick = 2.
						!p.charsize = 1.5
						!p.symsize = 2.
						!p.font = 2
						!x.title = 'Phase Difference'
						!y.title = 'Flare Maximum Flux'
						!x.range = [-0.5,0.5]
						!y.range = [0,MAX(flareflux)+0.05]
		
						PLOT, phasediff[WHERE(flareflux GE flarelowerlimit[a] AND flareflux LT flareupperlimit[a])], flareflux[WHERE(flareflux GE flarelowerlimit[a] AND flareflux LT flareupperlimit[a])], PSYM = 8

						CLEANPLOT, /SILENT
						DEVICE, /CLOSE
						SET_PLOT, 'x'
	
						IF N_ELEMENTS(flaremax) GT 1 AND N_ELEMENTS(phasediff[WHERE(flareflux GE flarelowerlimit[a] AND flareflux LT flareupperlimit[a] AND phasediff LT 10.)]) THEN BEGIN
							SET_PLOT, 'ps'
							DEVICE, FILENAME = directoryname + 'Flare_Compare_Hist/' + kepID[i]  + '_flare_flux' + flarelist[a] + '_jdfix' + jdlist[b] + '_endjdfix' + endjdlist[c] + '_2minmove_100bins_bw_histplot.eps', /LANDSCAPE, /COLOR, DECOMPOSED = 1
							!x.style = 1.
							!y.style = 1.
							!p.thick = 2.
							!x.thick = 2.
							!y.thick = 2.
							!p.charsize = 1.5
							!p.symsize = 2.
							!p.font = 2
							!x.title = 'Phase Difference'
							!y.title = 'Number of Flares'
			
							PLOTHIST, phasediff[WHERE(flareflux GE flarelowerlimit[a] AND flareflux LT flareupperlimit[a])], bin = 0.01, xrange = [-0.5,0.5]

							CLEANPLOT, /SILENT
							DEVICE, /CLOSE
							SET_PLOT, 'x'
						ENDIF
					ENDIF
				ENDFOR
			ENDIF
		ENDIF
	ENDFOR		

	FOR a = 0, N_ELEMENTS(flarelist) - 1 DO BEGIN
		toplefty = DBLARR(bins)	

		FOR d = 0, N_ELEMENTS(bottomleftx) - 1 DO BEGIN
			IF d NE N_ELEMENTS(bottomleftx) - 1 THEN BEGIN
				tempvals = (WHERE(allphasediffs GE -0.50+d/100. AND allphasediffs LT (-0.5+0.01)+d/100. AND allflarefluxes GE flarelowerlimit[a] AND allflarefluxes LT flareupperlimit[a]))
			ENDIF ELSE BEGIN ; the below is just for the end of the interval 
				tempvals = (WHERE(allphasediffs GE -0.50+d/100. AND allphasediffs LE (-0.5+0.01)+d/100. AND allflarefluxes GE flarelowerlimit[a] AND allflarefluxes LT flareupperlimit[a]))
			ENDELSE
			IF N_ELEMENTS(tempvals) GE 1 AND tempvals[0] GE 0. THEN BEGIN
				toplefty[d] = N_ELEMENTS(tempvals)
			ENDIF
		ENDFOR
		toprighty = toplefty

		SET_PLOT, 'ps'
		DEVICE, FILENAME = directoryname + 'All_stars_flare_flux' + flarelist[a] + '_jdfix' + jdlist[b] + '_endjdfix' + endjdlist[c] + '_2minmove_1spot_100bins_log_bw_plot.eps', /LANDSCAPE, /COLOR, DECOMPOSED = 1
		!x.style = 1.
		!y.style = 1.
		!p.thick = 2.
		!x.thick = 2.
		!y.thick = 2.
		!p.charsize = 1.5
		!p.symsize = 2.
		!p.font = 2	
		!x.title = textoidl('Phase Difference, \Delta\varphi')
		!y.title = 'Flare Maximum Flux Increase (Percent)'
		!x.range = [-0.5,0.5]
		!y.range = [0.95,115]
		IF N_ELEMENTS(WHERE(allflarefluxes GE flarelowerlimit[a] AND allflarefluxes LT flareupperlimit[a])) GT 1 THEN BEGIN
			PRINT, N_ELEMENTS(WHERE(allflarefluxes GE flarelowerlimit[a] AND allflarefluxes LT flareupperlimit[a]))
			PLOT, allphasediffs[WHERE(allflarefluxes GE flarelowerlimit[a] AND allflarefluxes LT flareupperlimit[a])], 100.*allflarefluxes[WHERE(allflarefluxes GE flarelowerlimit[a] AND allflarefluxes LT flareupperlimit[a])], PSYM = 8, /ylog
		ENDIF

		CLEANPLOT, /SILENT
		DEVICE, /CLOSE	
		SET_PLOT, 'x'

		SET_PLOT, 'ps'
		DEVICE, FILENAME = directoryname + 'All_stars_flare_flux' + flarelist[a] + '_jdfix' + jdlist[b] + '_endjdfix' + endjdlist[c] +  '_2minmove_1spot_100bins_bw_histplot.eps', /LANDSCAPE, /COLOR, DECOMPOSED = 1
		!x.style = 1.
		!y.style = 1.	
		!p.thick = 2.
		!x.thick = 2.
		!y.thick = 2.
		!p.charsize = 1.5
		!p.symsize = 2.
		!p.font = 2
		!x.title = textoidl('Phase Difference, \Delta\varphi')
		!y.title = 'Number of Flares'
		!x.range = [-0.5,0.5]
		!y.range = [0,MAX(toplefty)+5]

		PLOT, bottomleftx, bottomlefty, /NODATA
		FOR d = 0, N_ELEMENTS(bottomleftx) - 1 DO BEGIN
			polyfill, [bottomleftx[d],topleftx[d],toprightx[d],bottomrightx[d]],[bottomlefty[d],toplefty[d],toprighty[d],bottomrighty[d]], COLOR = 'CCCCCC'x	
		ENDFOR
		FOR d = 0, N_ELEMENTS(bottomleftx) - 1 DO BEGIN

			OPLOT, [bottomleftx[d],bottomrightx[d]], [bottomlefty[d],bottomrighty[d]], thick = 4, COLOR = '555555'x
			OPLOT, [bottomleftx[d],topleftx[d]], [bottomlefty[d],toplefty[d]], thick = 4, COLOR = '555555'x 
			OPLOT, [topleftx[d],toprightx[d]], [toplefty[d],toprighty[d]], thick = 4, COLOR = '555555'x
			OPLOT, [toprightx[d],bottomrightx[d]], [toprighty[d],bottomrighty[d]], thick = 4, COLOR = '555555'x

			IF a EQ 0 THEN BEGIN
				numflares[0,d] = MEAN([bottomleftx[d], bottomrightx[d]])
				numflares[1,d] = toplefty[d] 
			ENDIF 
			IF a EQ 1 THEN BEGIN
				numflares[2,d] = toplefty[d]
			ENDIF
		ENDFOR
					
		OPLOT, [-0.5,0.5], [0,0], thick = 3, COLOR = '000000'x
		OPLOT, [-0.5,0.5], [MAX(toplefty)+5,MAX(toplefty)+5], THICK = 3, COLOR = '000000'x
		OPLOT, [0.5,0.5], [0,MAX(toplefty)+5], THICK = 3, COLOR = '000000'x
		OPLOT, [-0.5,-0.5], [0,MAX(toplefty)+5], THICK = 3, COLOR = '000000'x



		CLEANPLOT, /SILENT
		DEVICE, /CLOSE
		SET_PLOT, 'x'

		PRINT, 'Actual starspot and flares count:  ', actualspotcount
		PRINT, 'Number of flares:  ', TOTAL(toplefty)

	ENDFOR




	SET_PLOT, 'ps'
	DEVICE, FILENAME = directoryname + 'All_stars_flare_flux' + flarelist[0] + '_' + flarelist[1] + '_jdfix' + jdlist[b] + '_endjdfix' + endjdlist[c] +  '_2minmove_1spot_100bins_color_histplot.eps', /LANDSCAPE, /COLOR, DECOMPOSED = 1
	!x.style = 1.
	!y.style = 1.	
	!p.thick = 2.
	!x.thick = 2.
	!y.thick = 2.
	!p.charsize = 1.5
	!p.symsize = 2.
	!p.font = 2
	!x.title = textoidl('Phase Difference, \Delta\varphi')
	!y.title = 'Number of Flares'
	!x.range = [-0.5,0.5]
	!y.range = [0,50]

	PLOT, bottomleftx, bottomlefty, /NODATA

	FOR a = 0, N_ELEMENTS(flarelist) - 1 DO BEGIN
		toplefty = DBLARR(bins)	

		FOR d = 0, N_ELEMENTS(bottomleftx) - 1 DO BEGIN
			tempvals = (WHERE(allphasediffs GE -0.50+d/100. AND allphasediffs LT (-0.5+0.01)+d/100. AND allflarefluxes GE flarelowerlimit[a] AND allflarefluxes LT flareupperlimit[a]))

			IF N_ELEMENTS(tempvals) GE 1 AND tempvals[0] GE 0. THEN BEGIN
				toplefty[d] = N_ELEMENTS(tempvals)
			ENDIF
		ENDFOR
		toprighty = toplefty



		IF a EQ 0 THEN BEGIN

		
			FOR d = 0, N_ELEMENTS(bottomleftx) - 1 DO BEGIN
				polyfill, [bottomleftx[d],topleftx[d],toprightx[d],bottomrightx[d]],[bottomlefty[d],toplefty[d],toprighty[d],bottomrighty[d]], COLOR = '89AFEF'x ;'CCCCCC'x	
			ENDFOR
			FOR d = 0, N_ELEMENTS(bottomleftx) - 1 DO BEGIN

				OPLOT, [bottomleftx[d],bottomrightx[d]], [bottomlefty[d],bottomrighty[d]], thick = 4, COLOR = '1360E2'x ;'555555'x
				OPLOT, [bottomleftx[d],topleftx[d]], [bottomlefty[d],toplefty[d]], thick = 4, COLOR = '1360E2'x ;'555555'x 
				OPLOT, [topleftx[d],toprightx[d]], [toplefty[d],toprighty[d]], thick = 4, COLOR = '1360E2'x ;'555555'x
				OPLOT, [toprightx[d],bottomrightx[d]], [toprighty[d],bottomrighty[d]], thick = 4, COLOR = '1360E2'x ;'555555'x
			ENDFOR


		ENDIF
		

		IF a EQ 1 THEN BEGIN

		
			FOR d = 0, N_ELEMENTS(bottomleftx) - 1 DO BEGIN
				polyfill, [bottomleftx[d],topleftx[d],toprightx[d],bottomrightx[d]],[bottomlefty[d],toplefty[d],toprighty[d],bottomrighty[d]], COLOR = 'AD83B8'x ;'CCCCCC'x	
			ENDFOR
			FOR d = 0, N_ELEMENTS(bottomleftx) - 1 DO BEGIN
	
				OPLOT, [bottomleftx[d],bottomrightx[d]], [bottomlefty[d],bottomrighty[d]], thick = 4, COLOR = '961CB6'x ;'555555'x
				OPLOT, [bottomleftx[d],topleftx[d]], [bottomlefty[d],toplefty[d]], thick = 4, COLOR = '961CB6'x ;'555555'x 
				OPLOT, [topleftx[d],toprightx[d]], [toplefty[d],toprighty[d]], thick = 4, COLOR = '961CB6'x ;'555555'x
				OPLOT, [toprightx[d],bottomrightx[d]], [toprighty[d],bottomrighty[d]], thick = 4, COLOR = '961CB6'x ;'555555'x
			ENDFOR


		ENDIF
					

	ENDFOR
		OPLOT, [-0.5,0.5], [0,0], thick = 3, COLOR = '000000'x
		OPLOT, [-0.5,0.5], [50,50], THICK = 3, COLOR = '000000'x
		OPLOT, [0.5,0.5], [0,50], THICK = 3, COLOR = '000000'x
		OPLOT, [-0.5,-0.5], [0,50], THICK = 3, COLOR = '000000'x
		CLEANPLOT, /SILENT
		DEVICE, /CLOSE
		SET_PLOT, 'x'


ENDFOR	
ENDFOR		


END