pro SAUCES_Kepler_missing_data, directoryname, target, gapsize
;
; This looks for the gaps in the data and records them
; so that those gaps can be used later on
; Also finds big jumps
; 
; R. Roettenbacher 8 June 2018
; 4 January 2019 Cleaned up for github

LC_CBV_filename = directoryname + 'FullLCs_CBV/' + target + '_LC_CBV.txt'

READCOL, LC_CBV_filename, lcdate, lcflux, FORMAT = 'd,d', /SILENT

OPENW, 1, directoryname + 'MissingData/' + target + '_MissingData.txt'

missing = 0

FOR i = 1, N_ELEMENTS(lcdate) - 1 DO BEGIN
	IF (lcdate[i] - lcdate[i-1]) GT gapsize THEN BEGIN
		PRINTF, 1, lcdate[i-1], lcdate[i]
		missing++
	ENDIF

ENDFOR

CLOSE, 1

RETURN 
END