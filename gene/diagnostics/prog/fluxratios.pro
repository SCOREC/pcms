; May 2018 tbg changes.
;   Major rewrite to ensure that diffusivities are
;   computed with actual, local gradients(!) and that
;   modern nrg files with up to 10 columns are properly
;   processed.
;   For simplicity, IV runs are assumed to have identical
;   last time stamps both in nrg and mom which is default
;   in GENE since 2017
;   A possible improvement would be to just store those
;   values that are computed. For now, the philosophy by
;   M. Kammerer has been kept/reintroduced.
;   WARNING: CHANGES IN THE SCAN PARAMETERS BY GENE ITSELF
;   ARE NOT TAKEN INTO ACCOUNT (CONCERNS KYMIN, FOR INSTANCE)

; Sep 9 2013 Ken Liao changes.
;   changed some help text
;   removed hardcoded limit on 70 entries in scan
;   changed phi to phi^2 in numerator and denominator
;   fix calculation of ratio_local2
;   fix ytitle and color for 2-dimensional surface plots
;   fix species label in comment line in ascii output
;   change ascii append file behavior
;   It is wrong to try to average over last 10 lines of nrg, since the
;   nrg file is organized such that each time point takes several
;   lines (depending on number of species). I have changed it to take
;   last data point, or the closest point to mom_time if need_phi =
;   true. This is only fixed for initial value problem solver.
FUNCTION fluxratios_info

  COMMON global_vars

  IF PTR_VALID(scanlog) THEN BEGIN
     ndimarr = N_ELEMENTS((*scanlog).dimarr)
     n_extvars = ndimarr+2
     ext_vars = STRARR(3,n_extvars)
     ext_vars[0,0] = 'vars'
     ext_vars[0,1] = 'species'
     ext_vars[0,2] = 'eigenvalues'
     ext_vars[0,3:n_extvars-1] = (*scanlog).dimarr[1:ndimarr-1]
     ext_vars[1,*] = '0'
     ext_vars[2,0] = 'array containing variable pairs, '+$
                     'e.g.: [4,0] gives <Gamma_es> / <|n|^2> '+$
                     '(default: [100,102])'
     ext_vars[2,1] = 'array containing pairs of '+$
                     'species for the variable pairs, where 0 is the '+$
                     'first species and -1 is replaced by all gui selected '+$
                     'species (default: [0,-1]'
     ext_vars[2,2:n_extvars-1] = 'index between '+rm0es(0)+' and '+$
        rm0es((*scanlog).dsize[*]-1)+' (corresponds to values from '+$
	rm0es(STRING([1,(*scanlog).axxis[*,0]]))+' to '+$
        rm0es(STRING([(*scanlog).dsize[0]-1,(*scanlog).rd[*]]))+'); default: * (all)'
  ENDIF ELSE ext_vars = ''

  RETURN, {$
          type      : 'scan',$
          title     : 'Flux ratio scan',$
          help_text : ['Plots the ratios between the selected variables. '+$
                       'One ratio is calculated per pair. '+$
                       'Variables 0 to 7 correspond to the nrg variables, '+$
                       '100 to D_es, 101 to D_em, 102 to chi_es, 103 to chi_em, '+$
                       '104 to phi^2, and 105 to phi*n'],$
          ext_vars  : ext_vars}

END

;######################################################################

PRO fluxratios_init, diag

  COMMON global_vars

  IF NOT PTR_VALID(scanlog) THEN BEGIN
     PRINT, 'skipping ' + (*diag).name + ': no scanfile found'
     RETURN
  ENDIF

  IF gui.out.n_spec_sel EQ 0 THEN BEGIN
     PRINT, 'No species selected!'
     RETURN
  ENDIF

  ; set default values for vars and species
  e = (*diag).external_vars
  n_extvars = N_ELEMENTS((*diag).table_entry[3,*])
  IF N_ELEMENTS(*e.vars) LT 1 THEN *e.vars = [100,102]
  IF N_ELEMENTS(*e.species) LT 1 THEN *e.species = [0,-1]
  ; now consider the scan range limitations (EV and scan dimensions)
  eigenvalues = (*diag).table_entry[3,2]
  IF STRLEN(eigenvalues) EQ 0 THEN eigenvalues='*' ELSE BEGIN
;     printerror, (*diag).name+': EIGENVALUES CURRENTLY RESTRICTED TO *'
;     eigenvalues='*'
  ENDELSE
  ;plotparA contains the scan ranges that will be considered
  ;for evaluation and plotting
  plotparA=REFORM((*diag).table_entry[3,3:n_extvars-1])
  inds = WHERE(plotparA EQ '')
  IF (inds[0] GE 0) THEN plotparA[inds]='*'

  ;SANITY CHECKS
  IF (N_ELEMENTS(*e.vars) MOD 2) THEN BEGIN
     printerror, (*diag).name+": Even numbers required for variable pairs"
     (*diag).selected = 0
     RETURN
  ENDIF

  variables = FIX(*e.vars)

  ; exception handling for false input:
  w = WHERE((variables LT 0) OR ((variables GT par.nrgcols) $
            AND (variables LT 100)) OR (variables GT 105))
  IF (w[0] NE -1) THEN BEGIN
     printerror, (*diag).name+': Incorrect input in vars'
     (*diag).selected = 0
     RETURN
  ENDIF

  ; exception handling for false species input:
  IF N_ELEMENTS(*e.species) NE N_ELEMENTS(*e.vars) THEN BEGIN
     printerror, (*diag).name+': number of species pairs does not match variable pairs'
     (*diag).selected = 0
     RETURN
  ENDIF
  species = FIX(*e.species)
  w = WHERE((species LT -1) OR (species GT gui.out.n_spec_sel))
  IF (w[0] NE -1) THEN BEGIN
     printerror, (*diag).name+': Incorrect input in spec'
     (*diag).selected = 0
     RETURN
  ENDIF

  ;now loop over species and repeat variable pairs if species = -1 is found
  first = 1
  FOR v = 0, N_ELEMENTS(*e.species)/2 - 1 DO BEGIN
     IF (*e.species)[2*v] GE 0 AND (*e.species)[2*v+1] GE 0 THEN BEGIN
        if (first) THEN BEGIN
           var_pairs = [variables[2*v:2*v+1]]
           spec_pairs = [species[2*v:2*v+1]]
           first = 0
        ENDIF ELSE BEGIN
           var_pairs = [[var_pairs],[variables[2*v:2*v+1]]]
           spec_pairs = [[spec_pairs],[species[2*v:2*v+1]]]
        ENDELSE
     ENDIF ELSE IF (*e.species)[2*v] EQ -1 AND (*e.species)[2*v+1] GE 0 THEN BEGIN
        FOR isp = 0, gui.out.n_spec_sel-1 DO BEGIN
           IF (first) THEN BEGIN
              var_pairs = [variables[2*v:2*v+1]]
              spec_pairs = [(*gui.out.spec_select)[isp],species[2*v+1]]
              first = 0
           ENDIF ELSE BEGIN
              var_pairs = [[var_pairs],[variables[2*v:2*v+1]]]
              spec_pairs = [[spec_pairs],[(*gui.out.spec_select)[isp],species[2*v+1]]]
           ENDELSE
        ENDFOR
     ENDIF ELSE IF (*e.species)[2*v] GE 0 AND (*e.species)[2*v+1] EQ -1 THEN BEGIN
        FOR isp = 0, gui.out.n_spec_sel-1 DO BEGIN
           IF (first) THEN BEGIN
              var_pairs = [variables[2*v:2*v+1]]
              spec_pairs = [species[2*v],(*gui.out.spec_select)[isp]]
              first = 0
           ENDIF ELSE BEGIN
              var_pairs = [[var_pairs],[variables[2*v:2*v+1]]]
              spec_pairs = [[spec_pairs],[species[2*v],(*gui.out.spec_select)[isp]]]
           ENDELSE
        ENDFOR
     ENDIF ELSE BEGIN
        FOR isp = 0, gui.out.n_spec_sel-1 DO BEGIN
           FOR isp2 = 0, gui.out.n_spec_sel-1 DO BEGIN
              IF (first) THEN BEGIN
                 var_pairs = [variables[2*v:2*v+1]]
                 spec_pairs = [(*gui.out.spec_select)[isp],(*gui.out.spec_select)[isp2]]
                 first = 0
              ENDIF ELSE BEGIN
                 var_pairs = [[var_pairs],[variables[2*v:2*v+1]]]
                 spec_pairs = [[spec_pairs],[(*gui.out.spec_select)[isp],(*gui.out.spec_select)[isp2]]]
              ENDELSE
           ENDFOR
        ENDFOR
     ENDELSE
  ENDFOR

  ;determine scan sizes and EV/IV mode
  ev_solv = par.n_ev GT 0  ; zero for initial value solver, 1 for eigenvalue solver
  nevs = (*scanlog).dsize[0] ;number of requested EVs (1 for IV)
  nscans = PRODUCT((*scanlog).dsize[1:*],/INTEGER)
  ndsize = PRODUCT((*scanlog).dsize,/INTEGER)
  nratios = N_ELEMENTS(var_pairs[0,*])

  w = WHERE((variables GE 104) AND (variables LE 105))
  need_phi = (w[0] NE -1)

  i = set_internal_vars(diag,{$
      var_pairs   : var_pairs,$
      spec_pairs  : spec_pairs,$
      eigenvalues : eigenvalues,$ ;exactly this name required for scan_lib.pro!!
      plotparA    : plotparA,$ ;contains scan ranges
      ev_solv     : ev_solv,$
      n_ev        : nevs,$
      need_phi    : need_phi,$
      phi_rms     : 0.0,$
      ratio_arr   : MAKE_ARRAY([nratios,nevs,nscans],/DOUBLE)})

  ;request nrg variables and - if needed - phi
  series.request_nrg = 1
  IF (need_phi) THEN fft_format, kxky=0

END

;######################################################################

FUNCTION fluxratios_get_var, diag, var, sp, tind

  COMMON global_vars

  i = (*diag).internal_vars

  IF (var LT par.nrgcols) THEN BEGIN ;nrg entries
     RETURN, (*nrg)[tind].data[var+sp*par.nrgcols]
  ENDIF ELSE BEGIN ;additional/special variables
     IF par.norm_flux_projection THEN BEGIN
        geo_fac = TOTAL(sqrt((*series.geom).gxx)*(*series.geom).jac_norm) / $
                  N_ELEMENTS((*series.geom).jac_norm)
     ENDIF ELSE BEGIN
        geo_fac = TOTAL((*series.geom).gxx*(*series.geom).jac_norm) / $
                  N_ELEMENTS((*series.geom).jac_norm)
     ENDELSE
     CASE var OF
        100 : BEGIN ; D_es
           RETURN, spec[sp].omn EQ 0.0 ? !VALUES.F_NAN : $
                   (*nrg)[tind].data[4+sp*par.nrgcols] * series.Lref / $
                   (spec[sp].omn * spec[sp].dens * series.nref * geo_fac)
        END
        101 : BEGIN ;D_em
           RETURN, spec[sp].omn EQ 0.0 ? !VALUES.F_NAN : $
                   (*nrg)[tind].data[5+sp*par.nrgcols] * series.Lref / $
                   (spec[sp].omn * spec[sp].dens * series.nref * geo_fac)
        END
        102 : BEGIN ;chi_es
           RETURN, spec[sp].omt EQ 0.0 ? !VALUES.F_NAN : $
                   (*nrg)[tind].data[6+sp*par.nrgcols] * series.Lref / $
                   (spec[sp].omt * spec[sp].dens * series.nref * $
                    spec[sp].temp * series.Tref * series.Qref * geo_fac)
        END
        103 : BEGIN ;chi_em
           RETURN, spec[sp].omt EQ 0.0 ? !VALUES.F_NAN : $
                   (*nrg)[tind].data[7+sp*par.nrgcols] * series.Lref / $
                   (spec[sp].omt * spec[sp].dens * series.nref * $
                    spec[sp].temp * series.Tref * series.Qref * geo_fac)

        END
        104 : BEGIN ;phi_rms^2
           RETURN, ((*i).phi_rms)^2
        END
        105 : BEGIN ;phi_rms*n_rms
           RETURN, (*i).phi_rms*SQRT((*nrg)[tind].data[0+sp*par.nrgcols])
        END
        ELSE : BEGIN
           printerror, 'invalid var in fluxratios_get_var'
           exit
        END
     ENDCASE
  ENDELSE
END
;######################################################################

FUNCTION fluxratios_get_var_name, var, sp, fancy=fancy

  COMMON global_vars

  sp_str = STRMID(spec[sp].name,0,1)
  IF KEYWORD_SET(fancy) THEN $
     sp_str = '!D'+sp_str+'!N' $
  ELSE sp_str = ','+sp_str

  IF (var LT par.nrgcols) THEN BEGIN ;nrg entries
     RETURN, get_nrg_string(var,fancy=fancy)+sp_str
  ENDIF ELSE BEGIN ;additional/special variables
     CASE var OF
        100 : RETURN, get_nrg_string(10,fancy=fancy)+sp_str ;D_es
        101 : RETURN, get_nrg_string(11,fancy=fancy)+sp_str ;D_em
        102 : RETURN, get_nrg_string(12,fancy=fancy)+sp_str ;chi_es
        103 : RETURN, get_nrg_string(13,fancy=fancy)+sp_str ;chi_em
        104 : RETURN, '!7U!6!U2!N' ;phi_rms^2
        105 : RETURN, '!7U!6!Drms!N * n!Drms!N'+sp_str ;phi_rms*n_rms
        ELSE : BEGIN
           printerror, 'invalid var in fluxratios_get_var_name'
           exit
        END
     ENDCASE
  ENDELSE
END

;######################################################################

PRO fluxratios_loop, diag

  COMMON global_vars

  i = (*diag).internal_vars
  iscan=FIX(series.run_labels) ;id of current scan file
  nevs = N_ELEMENTS((*i).ratio_arr[0,*,0])

; check whether value is supposed to be considered for
; further post-processing or neglected due to user defined
; scan parameter ranges
  calc=scan_plotnec(diag,[(*i).eigenvalues,(*i).plotparA],iscan)

  IF NOT calc THEN RETURN

  IF (*i).need_phi THEN BEGIN
     IF NOT (*i).ev_solv OR ((ABS(FIX(mom_time)) NE 0) AND $
                             (ABS(FIX(mom_time)) LE (*i).n_ev)) THEN BEGIN

                                ; averaging over kx and ky
        phi_rms_z = TOTAL((2.0*TOTAL(ABS($
                    (*mom[0,0].kxky)[*,1:par.nky0-1,*])^2,2)+ABS($
                    (*mom[0,0].kxky)[*,0,*])^2),1)
                                ; z averaging
        phi_rms_z = TEMPORARY(phi_rms_z) * (*series.geom).jac_norm
        (*i).phi_rms = SQRT(TOTAL(phi_rms_z)/par.nz0)
     ENDIF
  ENDIF

; calc nrgindex for run and set ratiofield
  IF par.comp_type EQ 'IV'  THEN BEGIN
     ;select last nrg time step
     tind = N_ELEMENTS((*nrg).time) - 1
     FOR v = 0, N_ELEMENTS((*i).var_pairs[0,*])-1 DO BEGIN
        (*i).ratio_arr[v,0,iscan-1] = $
           fluxratios_get_var(diag,(*i).var_pairs[0,v],(*i).spec_pairs[0,v],tind)/$
           fluxratios_get_var(diag,(*i).var_pairs[1,v],(*i).spec_pairs[1,v],tind)
     ENDFOR
  ENDIF ELSE BEGIN
     FOR isort = 0, (*scanlog).dsize[0]-1 DO BEGIN
        found = 0
        ind_sortf = (iscan-1)*(*scanlog).dsize[0]+isort
        ind_ev = (*scanlog).sortfield[ind_sortf]

        IF ((*nrg)[N_ELEMENTS(*nrg)-1].time GT 0.0) THEN $
           ind_ev = ABS(ind_ev)

        tind = WHERE((*nrg).time EQ ind_ev)

        IF (*i).need_phi THEN BEGIN
           IF (ind_ev EQ FIX(mom_time)) THEN found=ind_ev
        ENDIF ELSE found=(tind[0] NE -1)*ind_ev

        IF (found GT 0) THEN BEGIN
           FOR v = 0, N_ELEMENTS((*i).var_pairs[0,*])-1 DO BEGIN
              (*i).ratio_arr[v,found-1,iscan-1] = $
                 fluxratios_get_var(diag,(*i).var_pairs[0,v],(*i).spec_pairs[0,v],tind)/$
                 fluxratios_get_var(diag,(*i).var_pairs[1,v],(*i).spec_pairs[1,v],tind)
           ENDFOR
        ENDIF ELSE print,'stable or no EVs found'
     ENDFOR
  ENDELSE

END

;######################################################################

PRO fluxratios_output, diag

  COMMON global_vars

  i = (*diag).internal_vars
  npairs = N_ELEMENTS((*i).var_pairs[0,*])

  FOR v=0, npairs-1 DO BEGIN
     varname= fluxratios_get_var_name((*i).var_pairs[0,v],(*i).spec_pairs[0,v],/fancy)+'/'+$
              fluxratios_get_var_name((*i).var_pairs[1,v],(*i).spec_pairs[1,v],/fancy)

;     print, rm_idl_fmt(varname)
;     FOR iscan=0, N_ELEMENTS((*i).ratio_arr[0,0,*])-1 DO $
;        print, (*i).ratio_arr[v,*,iscan]
;     print, ''

     scan_plot, diag,varname,REFORM((*i).ratio_arr[v,*,*]),$
                'non', legend_pos=[0.01, 0.06, 0.99, 0.16], $
                pos=[0.1,0.26,0.99,0.9],$
                outset=(0+(v EQ 0)-(v EQ npairs-1)+2*(npairs EQ 1)),$
                multi=[0,1,1] ;,suffix='_pair_'+rm0es(v)
  ENDFOR


END
