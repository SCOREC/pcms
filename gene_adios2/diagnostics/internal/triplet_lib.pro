;#############################################################################
;#    this library contains almost all triplet-related internal functions    #
;#############################################################################

FUNCTION triplet_renorm, var, Lref=Lref
; The normalization should be as follows:
; (c_ref/L_ref) (n_ref^2 m_ref T_ref) (rho_ref/L_ref)^2
; (n_0j^2 B0 T_0j m_j) |g_j|^2/F_0j dv_|| d_mu
; summed over kx,ky and integrated over z.

; Note that n_0j, T_0j, m_j etc. are _NOT_ included in the energy arrays
; and -- technically -- still need to be considered (par.Xref)

  COMMON global_vars

  IF NOT KEYWORD_SET(Lref) THEN Lref = series.Lref

  rho_ref = SQRT(series.mref*series.Tref/series.Qref) / series.Bref
  c_ref = SQRT(series.Tref*series.Qref/series.mref)

  E_ref = series.nref^2 * series.mref * series.Qref * series.Tref * $
    (rho_ref/series.Lref)^2

  CASE var OF
    0	 : RETURN, E_ref * c_ref / series.Lref ; dE/dt
    ELSE : printerror, 'error in triplet_renorm (var ' + $
             rm0es(var) + ' not found)'
  ENDCASE

END

;##########################################################################

PRO destroy_triplet_struct

  COMMON global_vars

  IF PTR_VALID(triplet_time) THEN PTR_FREE, triplet_time
  IF PTR_VALID(triplet_data) THEN PTR_FREE, triplet_data

END

;##########################################################################

FUNCTION get_triplet_time, run, get_n_steps=get_n_steps, first=first, last=last

  COMMON global_vars

  file_triplet = get_file_path('triplet',run)

  IF NOT file_triplet.exist THEN RETURN, 0.0

  OPENR, triplet_lun, file_triplet.path, /GET_LUN, ERROR=err
  IF err NE 0 THEN printerror, 'strange error occured in get_triplet_time'
  IF EOF(triplet_lun) THEN BEGIN
    PRINT, 'empty triplet file'
    RETURN, 0.0
  ENDIF

  n_lines = FILE_LINES(file_triplet.path)
  IF n_lines LE 3 THEN BEGIN
    printerror, 'triplet file has wrong format'
    RETURN, 0.0
  ENDIF

  time = 0.0D

  IF KEYWORD_SET(first) THEN BEGIN
    SKIP_LUN, triplet_lun, 3, /LINES
    READF, triplet_lun, time
    FREE_LUN, triplet_lun
    time *= time_renorm(0)
    RETURN, time
  ENDIF

  IF KEYWORD_SET(last) THEN BEGIN
    SKIP_LUN, triplet_lun, n_lines - 1, /LINES
    READF, triplet_lun, time
    FREE_LUN, triplet_lun
    time *= time_renorm(0)
    RETURN, time
  ENDIF

  IF KEYWORD_SET(get_n_steps) THEN BEGIN
     n_steps = (n_lines - 3) / par.n_fields
     IF 10 * FIX(n_steps) NE FIX(10*n_steps) THEN PRINT, $
       'WARNING: triplet contains step number incompatible with n_fields'
     n_steps = FIX(n_steps)
     FREE_LUN, triplet_lun
     RETURN, n_steps
  ENDIF

END

;##########################################################################

PRO read_triplet, run, time=time, data=data, $
  get_steps=get_steps, n_steps=n_steps
; Reads the triplet_[run] file and makes the corresponding data available
; as a time array and data[triplet,pod,field,time].
; Alternatively, set /get_steps to return only the step number in n_steps.

  COMMON global_vars

  n_steps = 0

  IF N_ELEMENTS(run) NE 1 THEN run = 0

  triplet_file = get_file_path('triplet',run)

  exist = FILE_TEST(triplet_file.path,/READ,/REGULAR)
  IF NOT exist THEN BEGIN
    printerror, triplet_file.path + ' does not exist, skipping triplet'
    RETURN
  ENDIF

  n_lines = FILE_LINES(triplet_file.path)
  IF n_lines LE 3 THEN BEGIN
    printerror, triplet_file.path + ' has wrong format, skipping triplet'
    RETURN
  ENDIF

  OPENR, triplet_lun, triplet_file.path, /GET_LUN, ERROR=err
  IF err NE 0 THEN BEGIN
    printerror, triplet_file.path + ' could not be read, skipping triplet'
    RETURN
  ENDIF

  header = STRARR(3)
  READF, triplet_lun, header

  n_steps = (n_lines - 3.0) / par.n_fields
  IF 10 * LONG(n_steps) NE LONG(10*n_steps) THEN BEGIN
    printerror, triplet_file.path + ' has wrong structure, skipping triplet'
    RETURN
  ENDIF ELSE n_steps = LONG(n_steps)
  IF KEYWORD_SET(get_steps) THEN BEGIN
    FREE_LUN, triplet_lun
    RETURN
  ENDIF

  rawdata = FLTARR(1+1L*par.n_triplets*par.n_pod_triplets,$
    n_steps*par.n_fields,/NOZERO)
  READF, triplet_lun, rawdata

  FREE_LUN, triplet_lun

  time = REFORM(rawdata[0,par.n_fields*LINDGEN(n_steps)])
  data = REFORM(rawdata[1:*,*],[par.n_pod_triplets,par.n_triplets,$
    par.n_fields,n_steps])

END

;##########################################################################

PRO get_triplet_data

  COMMON global_vars

  n_steps_all = LONARR(series.n_runs)
  FOR run = 0, series.n_runs - 1 DO BEGIN
    read_triplet, run, /get_steps, n_steps=n_steps
    n_steps_all[run] = n_steps
  ENDFOR

  n_steps_tot = TOTAL(n_steps_all)
  IF n_steps_tot LT 1 THEN RETURN
  start_step = series.n_runs GT 1 ? [0L,n_steps_all[0:series.n_runs-2]] : 0L

  triplet_time = PTR_NEW(pf_arr([n_steps_tot]),/NO_COPY)
  triplet_data = $
    PTR_NEW(pf_arr([par.n_pod_triplets,par.n_triplets,$
    par.n_fields,n_steps]),/NO_COPY)

  FOR run = 0, series.n_runs - 1 DO BEGIN
    read_triplet, run, time=time, data=data

    (*triplet_time)[start_step[run]] = time
    (*triplet_data)[0,0,0,start_step[run]] = data
  ENDFOR

END

;##########################################################################

PRO triplet_loop, diag0

  COMMON global_vars

  series.step_count = 1
  get_triplet_data

END
