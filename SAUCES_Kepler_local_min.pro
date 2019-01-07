pro SAUCES_Kepler_local_min, directoryname, target, period, FLATWRM, typeofstar
; 
;  This looks for the local min in a light curve
;  The goal will be to compare this to flare timing
;
;  The inputs are the target name, the rotation period, and whether or not 
;  there needs to be special output for the command to run FLATWRM
;  Use 'Y' to make the special output. 
;
;  The output are
;  --list of minima detected printed with the CBV-removed flux, and a smoothed flux
;
;
;  R. Roettenbacher 22 August 2016
;  Updated 4 January 2019 to be uploaded to github

LC_CBV_filename = directoryname + 'FullLCs_CBV/' + target + '_LC_CBV.txt'

READCOL, LC_CBV_filename, CBV_mbjd, CBV_flux, FORMAT = 'd,d'

PLOT, CBV_mbjd, CBV_flux, /NODATA
PLOT, CBV_mbjd, CBV_flux, XRANGE = [MIN(CBV_mbjd), MAX(CBV_mbjd)], YRANGE = [MIN(CBV_flux), MAX(CBV_flux)], XSTY = 1
PLOT, CBV_mbjd, CBV_flux, XRANGE = [1540, 1600], YRANGE = [MIN(CBV_flux), MAX(CBV_flux)], XSTY = 1
PLOT, CBV_mbjd, CBV_flux, XRANGE = [540, 600], YRANGE = [MIN(CBV_flux), MAX(CBV_flux)], XSTY = 1


min_time_filename = directoryname + 'TimeofMinima_CBV/' + target + '_TimeofMinima_CBV.txt' 

orig_and_smooth_lc_filename = directoryname + 'FullLCs_CBV_AND_SMOOTH/' + target + '_FullLCs_CBV_AND_SMOOTH.txt' 


FLATWRM_filename = directoryname + 'FLATWRM_input/' + target + '_FLATWRM_input.txt'


IF period GT 0.5 THEN BEGIN
	PRINT, 'Period above 0.5 days'

	fwhmcoeff = 0.1
	smooth_flux = most_gsmooth(CBV_flux, FWHM = fwhmcoeff*period, TIME = CBV_mbjd)

	PLOT, CBV_mbjd, smooth_flux, XRANGE = [MIN(CBV_mbjd), MAX(CBV_mbjd)], YRANGE = [MIN(smooth_flux), MAX(smooth_flux)], XSTY = 1
	PLOT, CBV_mbjd, smooth_flux, XRANGE = [1540,1600], YRANGE = [MIN(smooth_flux), MAX(smooth_flux)],  XSTY = 1
	oplot, CBV_mbjd, CBV_flux, linestyle = 1

	OPENW, 1, min_time_filename

	PRINTF, 1, '; Times of minima in Kepler light curves for ' + target + '.'
	PRINTF, 1, '; Period = ' + string(period) + ' days.'
	PRINTF, 1, '; Used CBV-removed light curve.'
	PRINTF, 1, '; Column 1:  Time of minimum as determined from the'
	PRINTF, 1, ';            smoothed light curve. '
	PRINTF, 1, '; Column 2:  Flux at the time of minimum from the'
	PRINTF, 1, ';            smoothed light curve. '
	PRINTF, 1, '; Column 3:  Flux at the time of minimum from the'
	PRINTF, 1, ';            ORIGINAL, UNSMOOTHED light curve.'

	FOR i = 10, N_ELEMENTS(CBV_mbjd) - 11 DO BEGIN
	
		; Check to see if the preceding and following data points are higher in flux 
		IF MIN(smooth_flux[i-10:i+10]) EQ smooth_flux[i] THEN BEGIN
			tempx = [CBV_mbjd[i], CBV_mbjd[i]]
			tempy = [0.,2.]
			OPLOT, tempx, tempy		
	
			PRINTF, 1, CBV_mbjd[i], smooth_flux[i], CBV_flux[i]

		ENDIF	

	ENDFOR

	OPENW, 3, orig_and_smooth_lc_filename

	PRINTF, 3, '; CBV-applied AND smoothed Kepler light curves for ' + target + '.'
	PRINTF, 3, '; Period = ' + string(period) + ' days.'
	PRINTF, 3, '; Column 1:  MBJD'
	PRINTF, 3, '; Column 2:  Flux, CBV-applied, SMOOTHED'
	PRINTF, 3, ';            light curve. '
	PRINTF, 3, '; Column 3:  Flux, ORIGINAL, CBV-applied, UNSMOOTHED'
	PRINTF, 3, ';            light curve.'
	
	FOR j = 0, N_ELEMENTS(CBV_mbjd) - 1 DO BEGIN
		PRINTF, 3, CBV_mbjd[j], smooth_flux[j], CBV_flux[j]
	ENDFOR
	CLOSE, 3


ENDIF ELSE BEGIN
	PRINT, 'Period below 0.5 days'

        fwhmcoeff = 0.2
        smooth_flux = most_gsmooth(CBV_flux, FWHM = fwhmcoeff*period, TIME = CBV_mbjd)

        PLOT, CBV_mbjd, smooth_flux, XRANGE = [MIN(CBV_mbjd), MAX(CBV_mbjd)], YRANGE = [MIN(smooth_flux), MAX(smooth_flux)], XSTY = 1
        PLOT, CBV_mbjd, smooth_flux, XRANGE = [1540,1600], YRANGE = [MIN(smooth_flux), MAX(smooth_flux)],  XSTY = 1
	oplot, CBV_mbjd, CBV_flux, linestyle = 1

        OPENW, 1, min_time_filename

	PRINTF, 1, '; Times of minima in Kepler light curves for ' + target + '.'
	PRINTF, 1, '; Period = ' + string(period) + ' days.'
	PRINTF, 1, '; Used CBV-removed light curve.'
	PRINTF, 1, '; Column 1:  Time of minimum as determined from the'
	PRINTF, 1, ';            smoothed light curve. '
	PRINTF, 1, '; Column 2:  Flux at the time of minimum from the'
	PRINTF, 1, ';            smoothed light curve. '
	PRINTF, 1, '; Column 3:  Flux at the time of minimum from the'
	PRINTF, 1, ';            ORIGINAL, UNSMOOTHED light curve.'

        FOR i = 5, N_ELEMENTS(CBV_mbjd) - 6 DO BEGIN

                ; Check to see if the preceding and following data points are higher in flux 
                IF MIN(smooth_flux[i-5:i+5]) EQ smooth_flux[i] THEN BEGIN
                        tempx = [CBV_mbjd[i], CBV_mbjd[i]]
                        tempy = [0.,2.]
                        OPLOT, tempx, tempy

                        PRINTF, 1, CBV_mbjd[i], smooth_flux[i], CBV_flux[i]

                ENDIF

        ENDFOR

	OPENW, 3, orig_and_smooth_lc_filename

	PRINTF, 3, '; CBV-applied AND smoothed Kepler light curves for ' + target + '.'
	PRINTF, 3, '; Period = ' + string(period) + ' days.'
	PRINTF, 3, '; Column 1:  MBJD'
	PRINTF, 3, '; Column 2:  Flux, CBV-applied, SMOOTHED'
	PRINTF, 3, ';            light curve. '
	PRINTF, 3, '; Column 3:  Flux, ORIGINAL, CBV-applied, UNSMOOTHED'
	PRINTF, 3, ';            light curve.'
	
	FOR j = 0, N_ELEMENTS(CBV_mbjd) - 1 DO BEGIN
		PRINTF, 3, CBV_mbjd[j], smooth_flux[j], CBV_flux[j]
	ENDFOR
	CLOSE, 3


ENDELSE

CLOSE, 1


IF FLATWRM EQ 'Y' OR FLATWRM EQ 'y' THEN BEGIN
	OPENW, 2, FLATWRM_filename
	PRINTF, 2, 'python ./ml_flare.py --flarepoints=2 -o ' + directoryname + 'FLATWRM_flarepts2_output/' + target + '_flarepts2_output.txt ' + directoryname + 'FullLCs_CBV_noheader/' + target + '_LC_CBV.txt'
	PRINTF, 2, 'python ./ml_flare.py --flarepoints=3 -o ' + directoryname + 'FLATWRM_flarepts3_output/' + target + '_flarepts3_output.txt ' + directoryname + 'FullLCs_CBV_noheader/' + target + '_LC_CBV.txt'
	CLOSE, 2
ENDIF

RETURN
END




