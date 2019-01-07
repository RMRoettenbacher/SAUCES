pro SAUCES_Kepler_CBV_remove, directoryname, keplerarchive, target
;
; This program is going to open and read in kepler fits files
; by quarter, which then have the CBVs removed 
; The output are files containing the entire Kepler light curve.
; --list of CBVs used for each star
; --full Kepler light curve, median-divided and CBVs removed
; --full Kepler light curve, median-divided and CBVs NOT removed
;
; R. Roettenbacher 8 Aug 2016
;
; Handling of CBVs adapted from B. Jackson's kepcotrend.pro, which is based off of the
; python kepcotrend program from PyKE.
; 
; Updated 3 January 2019 Cleaned and commented for github
; 
; This procedure requires being passed a directory name, the directory of the Kepler archive, 
; and a target (name as a string).


CLOSE, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10


; Read in the list of Kepler quarter start dates.  These are used to name
; each quarter of the data.  This file can also be found on github, if needed.
READCOL, directoryname + 'Kepler_quarter_suffixes.txt', quarterdatelist, FORMAT = 'a'

; Read in the list of cotrending basis vectors -- Note Q4 is listed twice because of two different start dates
; The ones included in this file are from data release 25.  This file can be found on github.  The 
; CBVs can be found https://archive.stsci.edu/kepler/cbv.html
READCOL, directoryname + 'Kepler_CBV_list.txt', CBVlist, FORMAT = 'a'

; Rewrite all of the cbvlist entries to include the whole directory name
; Edit as necessary for where you have your CBV files
CBVlist = directoryname + CBVlist

CBVsused = INTARR(19) ; 18 is the number of quarters -- 0-17 -- +1 extra for Q4 being repeated
CBVschisq = DBLARR(19)

; Name the output files
LC_CBV_filename = directoryname + 'FullLCs_CBV/' + target + '_LC_CBV.txt'
LC_CBV_filename_noheader = directoryname + 'FullLCs_CBV_noheader/' + target + '_LC_CBV_noheader.txt'
LC_PDC_filename = directoryname + 'FullLCs_PDC/' + target + '_LC_PDC.txt'
CBVs_used_filename = directoryname + 'CBVsUsed/' + target + '_CBVs.txt'

FILE_MKDIR, directoryname + 'QuarterLCs_CBV/' + target
FILE_MKDIR, directoryname + 'QuarterLCs_PDC/' + target


; Need a counter to make sure the light curve gets written correctly
k = 0

; To get the right light curves, need to designate the appropriate folders
bigfolder = strmid(target, 0, 7) 	; This is the first folder in the Kep Archive
smallfolder = strmid(target, 0, 8)	; This is the second folder in the Kep Archive
; This assumes a structure of 'kplr000/kplr0000/Q00/' for your Kepler archive

; Because Q4 has two dates, need to just hard code this as an array
quarternum = ['00','01','02','03','04','04','05','06','07','08','09','10','11','12','13','14','15','16','17']


; Locate the .fits files in the archive (target passed from target list in script)
; Go through each quarter --- note that Q4 has two values.  
FOR j = 0, 18 DO BEGIN
	cd, keplerarchive + bigfolder + '/' + smallfolder + '/Q' + quarternum[j] + '/'

	LC_CBV_filename_Q = directoryname + 'QuarterLCs_CBV/' + target + '/' + target + '_LC_CBV_Q' + quarternum[j] + '.txt'
	LC_PDC_filename_Q = directoryname + 'QuarterLCs_PDC/' + target + '/' + target + '_LC_PDC_Q' + quarternum[j] + '.txt'

	lcname =  target + quarterdatelist[j]

	; Skip the quarters that don't have data
	IF FILE_TEST(lcname) EQ 1 THEN BEGIN
			
		; If the quarter has data for a particular star, 
		; store info from the structure and header
		LCStruct = MRDFITS(lcname, 1, header, /SILENT)
		mbjd = DOUBLE(LCStruct.TIME)
		sapflux = DOUBLE(LCStruct.SAP_FLUX)			
		PDCflux = DOUBLE(LCStruct.PDCSAP_FLUX)
		sapqual = DOUBLE(LCStruct.SAP_QUALITY)

		LCRead = READFITS(lcname, readheader, /SILENT)
		modnum = SXPAR(readheader, 'MODULE')
		outnum = SXPAR(readheader, 'OUTPUT')

		; Location of the star on the Kepler chip...
		CBVmodout = 'MODOUT_' + STRTRIM(modnum,1) + '_' + STRTRIM(outnum,1) 	

		; Pull out the structure that contains the CBVs for the particular
		; module and output of the star 
		CBVStruct = MRDFITS(CBVlist[j], CBVmodout, header, /SILENT)

		CBVs = [[CBVStruct.VECTOR_1], [CBVStruct.VECTOR_2], [CBVStruct.VECTOR_3], [CBVStruct.VECTOR_4], [CBVStruct.VECTOR_5], [CBVStruct.VECTOR_6], [CBVStruct.VECTOR_7], [CBVStruct.VECTOR_8], [CBVStruct.VECTOR_9], [CBVStruct.VECTOR_10], [CBVStruct.VECTOR_11], [CBVStruct.VECTOR_12], [CBVStruct.VECTOR_13], [CBVStruct.VECTOR_14], [CBVStruct.VECTOR_15], [CBVStruct.VECTOR_16]]
		CBVin = INDGEN(16)  ; There are 16 CBVs for each quarter

		; There are a lot of NaNs in the Kepler light curves. This removes them (first and third parts).
		; There are a lot of flags for things like cosmic rays and Auger brightening.  This removes them (second part).	
		goodval = WHERE(sapflux GT 0. AND sapqual EQ 0 AND PDCflux GT 0.)

	        ; Get CBVs in a format agreeable with the form needed for the kepcotrend 
		; AND cut out the bad values
		CBVs = CBVs[goodval,*]
                CBVs = DOUBLE(CBVs[*,CBVin])

		; Update the sapflux to be lacking those NaNs
		sapflux = sapflux[goodval]
		mbjd = mbjd[goodval]	

		PDCflux = PDCflux[goodval]


		; Prepare for CBV removal by median subtraction and median division
		; Source:  To use the CBVs for least-squares fitting, subtract the 
		; median uncorrected flux from the uncorrected flux time series of 
		; interest and divide by the median. Since the basis is orthonormal, 
		; the linear least-squares fit coefficient of the nth CBV is simply 
		; the inner product of the median-removed, median-normalized uncorrected 
		; flux time series with the nth CBV. Subtract the fit to get the corrected 
		; (median-removed, median-normalized) flux time series.
		; Kepler Data Characteristics Handbook (17 Aug 2011) 
 		; http://archive.stsci.edu/kepler/manuals/Data_Characteristics_Handbook_20110817.pdf

		medadjsapflux = (sapflux - MEDIAN(sapflux))/MEDIAN(sapflux)

		; You can silence the plotting if you want.  
		!P.MULTI = [0,1,3]

		; Formula take form kepcotrend.py
		CBVcoeffs = MATRIX_MULTIPLY(INVERT(MATRIX_MULTIPLY(CBVs, CBVs, /ATRANSPOSE)), MATRIX_MULTIPLY(CBVs, medadjsapflux, /ATRANSPOSE))

		; Apply the CBVs, one at a time (adding on another with each step)
		FOR l = 0, 15 DO BEGIN 

			IF l EQ 0 THEN BEGIN
				previousflux = medadjsapflux
				CBVflux = medadjsapflux - CBVcoeffs[l]*CBVs[*,l]
			ENDIF ELSE BEGIN
				previousflux = CBVflux
				CBVflux = CBVflux - CBVcoeffs[l]*CBVs[*,l]
			ENDELSE

			chisq = TOTAL(((CBVflux - previousflux)^2.)/previousflux)

			; See what slope a line fit through the data would give.
			; This helps make sure that the removal removed enough CBVs
			linecoeffs = LINFIT(mbjd, CBVflux)
			slope = linecoeffs[1]
			intercept = linecoeffs[0]

			; Select the lowest number of CBVs that is only a small improvement
			; to the next CBV
			; Want to store the lowest CBV, so not the one just removed, but one
			; less.  In this if statement, identify that CBV and then revert
			; the light curve to it.

			IF l GT 0 AND ABS(chisq) LT 0.3 AND 1./ABS(chisq) GT 0. AND ABS(intercept) LT 1. THEN BEGIN ;AND ABS(slope) LT 2.d-4 AND THEN BEGIN  ; This criterion might also be necessary

				CBVsused[j] = l ; CBV desired isn't the one just used, but -1 
						; CBV number is offset by one, so l is the 
						; right number
				CBVschisq[j] = chisq

				; Need to add back the last CBV!
				CVflux = CBVflux + CBVcoeffs[l]*CBVs[*,l]

				; Undo the normalization and shifting
				; To go back to the original SAP flux "level" in order to 
				; remove the offset --- to the median divide below! 
				CBVflux = (CBVflux*MEDIAN(sapflux)) + MEDIAN(sapflux)
				k++ 		; counter for writing light curve
				l = 16		; kick out of loop
	

			ENDIF
		ENDFOR

		; Each quarter is individually normalized.

		IF CBVsused[j] NE 0 THEN BEGIN
			IF k EQ 1 THEN BEGIN
				FullLC_mbjd = mbjd
				FullLC_flux = CBVflux/MEDIAN(CBVflux)
				FullLC_pdc = PDCflux/MEDIAN(PDCflux)
			ENDIF ELSE BEGIN	
				FullLC_mbjd = [FullLC_mbjd,mbjd]
				FullLC_flux = [FullLC_flux,CBVflux/MEDIAN(CBVflux)]
				FullLC_pdc = [FullLC_PDC,PDCflux/MEDIAN(PDCflux)]
			ENDELSE
		ENDIF 


		OPENW, 3, LC_CBV_filename_Q
		PRINTF, 3, '; Kepler light curve for ' + target + '.'
		PRINTF, 3, '; Quarter ' + quarternum[j]
		PRINTF, 3, '; CBVs have been removed.'
		PRINTF, 3, '; Column 1:  Time from the Kepler .fits files (MBJD)'
		PRINTF, 3, '; Column 2:  Normalized flux with CBVs removed'
		FOR p = 0, N_ELEMENTS(mbjd) - 1 DO BEGIN
			PRINTF, 3, mbjd[p], FORMAT = '(F15.7, $)'
			PRINTF, 3, CBVflux[p]
		ENDFOR
		CLOSE, 3

		OPENW, 4, LC_PDC_filename_Q
		PRINTF, 4, '; Kepler light curve for ' + target + '.'
		PRINTF, 4, '; Quarter ' + quarternum[j]
		PRINTF, 4, '; CBVs have NOT been removed.'
		PRINTF, 4, '; Column 1:  Time from the Kepler .fits files (MBJD)'
		PRINTF, 4, '; Column 2:  Normalized flux without CBVs removed'
		FOR p = 0, N_ELEMENTS(mbjd) - 1 DO BEGIN
			PRINTF, 4, mbjd[p], FORMAT = '(F15.7, $)'
			PRINTF, 4, PDCflux[p]
		ENDFOR
		CLOSE, 4
	ENDIF
ENDFOR


PLOT, FullLC_mbjd, FullLC_flux, XRANGE = [MIN(FullLC_mbjd), MAX(FullLC_mbjd)],  YRANGE = [MIN(FullLC_flux), MAX(FullLC_flux)], XSTY = 1

OPENW, 1, CBVs_used_filename
PRINTF, 1, '; CBVs used on ' + target + '.'
PRINTF, 1, '; Note that the fifth entry is for the Q4 that uses the '
PRINTF, 1, '; anomalous start date.'
PRINTF, 1, '; Column 1:  Quarter number (first instance of Q4 is anomaly)'
PRINTF, 1, '; Column 2:  Number of CBVs used for that quarter light curve'
PRINTF, 1, '; Column 3:  chi-square'
FOR p = 0, 18 DO BEGIN
	PRINTF, 1, quarternum[p], CBVsused[p], CBVschisq[p]
ENDFOR
CLOSE, 1



OPENW, 1, LC_CBV_filename
OPENW, 2, LC_PDC_filename
OPENW, 5, LC_CBV_filename_noheader

PRINTF, 1, '; Kepler light curve for ' + target + '.'
PRINTF, 1, '; CBVs have been removed.'
PRINTF, 1, '; Column 1:  Time from the Kepler .fits files (MBJD)'
PRINTF, 1, '; Column 2:  Normalized flux with CBVs removed'

PRINTF, 2, '; Kepler light curve for ' + target + '.'
PRINTF, 2, '; CBVs have NOT been removed.'
PRINTF, 2, '; Column 1:  Time from the Kepler .fits files (MBJD)'
PRINTF, 2, '; Column 2:  Normalized PDC flux without CBVs removed'  


FOR m = 0, N_ELEMENTS(FullLC_mbjd) - 1 DO BEGIN
	PRINTF, 1, FullLC_mbjd[m], Format='(F15.7, $)'
	PRINTF, 1, FullLC_flux[m]
	PRINTF, 2, FullLC_mbjd[m], Format='(F15.7, $)'
	PRINTF, 2, FullLC_pdc[m]
	PRINTF, 5, FullLC_mbjd[m], Format='(F15.7, $)'
	PRINTF, 5, FullLC_flux[m]
ENDFOR
CLOSE, 1
CLOSE, 2	

RETURN
END
