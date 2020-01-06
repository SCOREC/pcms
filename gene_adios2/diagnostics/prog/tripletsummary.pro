FUNCTION tripletsummary_info

  RETURN, {$
    type      : 'triplet',$
    title     : 'Triplet summary',$
    help_text : ['Plots time traces of fields from the triplet file.'],$
    ext_vars  : [['pod','0','POD number(s); default: all'],$
                 ['triplet_ind','0','These triplets '+$
                  'will be plotted; default: all'],$
                 ['field','0','Fields to be plotted (Phi,A_par,B_par); '+$
                  'default: all'],$
                 ['group_plots','0','How to group data into plots: '+$
                  '0: all dimensions separately; 1: combine all POD; '+$
                  '2: combine POD, triplets (default); 3: combine all'],$
                 ['list_triplets','1','List all available triplets']]}

END

;######################################################################

PRO tripletsummary_init, diag

  COMMON global_vars

  e = (*diag).external_vars
  pod = N_ELEMENTS(*e.pod) LT 1 ? INDGEN(par.n_pod_triplets) : *e.pod
  triplet_ind = N_ELEMENTS(*e.triplet_ind) LT 1 ? INDGEN(par.n_triplets) : $
    *e.triplet_ind
  field = N_ELEMENTS(*e.field) LT 1 ? INDGEN(par.n_fields) : *e.field
  group_plots = N_ELEMENTS(*e.group_plots) NE 1 ? 2 : *e.group_plots
  list_triplets = KEYWORD_SET(*e.list_triplets)

  i = set_internal_vars(diag,$
    {pod           : pod,$
     triplet_ind   : triplet_ind,$
     field         : field,$
     group_plots   : group_plots,$
     list_triplets : list_triplets})

END

;######################################################################

PRO tripletsummary_loop, diag

END

;######################################################################

PRO tripletsummary_output, diag

  COMMON global_vars

  i = (*diag).internal_vars

  IF N_ELEMENTS(*triplet_time) LE 1 THEN BEGIN
    PRINT, (*diag).name + ': less than two time steps, skipping output'
    RETURN
  ENDIF

  tinds = WHERE((*triplet_time GE gui.out.start_t) AND $
    (*triplet_time LE gui.out.end_t),tcount)
  IF tcount LT 2 THEN BEGIN
    PRINT, (*diag).name + ': < 2 time steps, skipping output'
    RETURN
  ENDIF

  n_pod = N_ELEMENTS((*i).pod)
  n_triplets = N_ELEMENTS((*i).triplet_ind)
  n_fields = N_ELEMENTS((*i).field)

  time = (*triplet_time)[tinds]
  uniq_inds = UNIQ(time)
  time = time[uniq_inds]

  triplet_dat = (*triplet_data)[*,*,*,tinds]
  triplet_dat = triplet_dat[*,*,(*i).field,*]
  triplet_dat = triplet_dat[*,(*i).triplet_ind,*,*]
  triplet_dat = triplet_dat[(*i).pod,*,*,*]
  triplet_dat = REFORM(triplet_dat,[n_pod,n_triplets,n_fields,tcount],$
    /OVERWRITE)
  tcount = N_ELEMENTS(uniq_inds)
  triplet_dat = REFORM(triplet_dat[*,*,*,uniq_inds],$
    [n_pod,n_triplets,n_fields,tcount])

  set_output, diag, /ps

  CASE (*i).group_plots OF
    0 : BEGIN
      FOR f = 0, n_fields - 1 DO FOR t = 0, n_triplets - 1 DO $
        FOR p = 0, n_pod - 1 DO PLOT, time, triplet_dat[p,t,f,*], COLOR=1, $
        /XSTYLE, /YSTYLE, XTITLE=get_var_string(0,/time,/units,/fancy), $
        YTITLE='!6amplitude', TITLE='!6'
    END
    1 : BEGIN
      FOR f = 0, n_fields - 1 DO FOR t = 0, n_triplets - 1 DO BEGIN
        maxval = MAX(triplet_dat[*,t,f,*],MIN=minval)
        PLOT, [time[0],time[tcount-1]], [minval,maxval], COLOR=1, /NODATA, $
          /XSTYLE, /YSTYLE, XTITLE=get_var_string(0,/time,/units,/fancy), $
          YTITLE='!6amplitude', TITLE='!6'

        FOR p = 0, n_pod - 1 DO OPLOT, time, triplet_dat[p,t,f,*], COLOR=p
      ENDFOR
    END
    2 : BEGIN
      FOR f = 0, n_fields - 1 DO BEGIN 
        maxval = MAX(triplet_dat[*,*,f,*],MIN=minval)
        PLOT, [time[0],time[tcount-1]], [minval,maxval], COLOR=1, /NODATA, $
          /XSTYLE, /YSTYLE, XTITLE=get_var_string(0,/time,/units,/fancy), $
          YTITLE='!6amplitude', TITLE='!6'

        FOR t = 0, n_triplets - 1 DO FOR p = 0, n_pod - 1 DO $
          OPLOT, time, triplet_dat[p,t,f,*], COLOR=p, LINE=t
      ENDFOR
    END
    3 : BEGIN
      maxval = MAX(triplet_dat,MIN=minval)
      PLOT, [time[0],time[tcount-1]], [minval,maxval], COLOR=1, /NODATA, $
        /XSTYLE, /YSTYLE, XTITLE=get_var_string(0,/time,/units,/fancy), $
        YTITLE='!6amplitude', TITLE='!6'

      FOR f = 0, n_fields - 1 DO FOR t = 0, n_triplets - 1 DO $
        FOR p = 0, n_pod - 1 DO OPLOT, time, triplet_dat[p,t,f,*], $
        COLOR=(p+t*n_pod<255), LINE=f
    END
    ELSE :
  ENDCASE

  IF (*i).list_triplets THEN BEGIN
    n_triplets_pp = 14
    n_pages = CEIL(1.0*par.n_pod_triplets/n_triplets_pp)

    xcoord = [0.1,0.3,0.45,0.65,0.8]
    header = ['!6i!Dtriplet!N','!6k!S!Dx!N!R!U(1)!N',$
      '!6k!S!Dy!N!R!U(1)!N','!6k!S!Dx!N!R!U(2)!N','!6k!S!Dy!N!R!U(2)!N']

    kx_triplets = (*series.kx)[((*par.kx_triplets)[3*INDGEN(par.n_triplets)] + $
      par.nx0) MOD par.nx0,0]
    ky_triplets = (*series.ky)[(*par.ky_triplets)[3*INDGEN(par.n_triplets)]]
    kxp_triplets = (*series.kx)[((*par.kx_triplets)[3*INDGEN(par.n_triplets)+1] + $
      par.nx0) MOD par.nx0,0]
    kyp_triplets = (*series.ky)[(*par.ky_triplets)[3*INDGEN(par.n_triplets)+1]]

    triplets_str = STRARR(par.n_triplets,5)
    triplets_str[0,0] = rm0es(INDGEN(par.n_triplets))
    triplets_str[0,1] = rm0es(kx_triplets,prec=4)
    triplets_str[0,2] = rm0es(ky_triplets,prec=4)
    triplets_str[0,3] = rm0es(kxp_triplets,prec=4)
    triplets_str[0,4] = rm0es(kyp_triplets,prec=4)

    FOR p = 0, n_pages - 1 DO BEGIN
      PLOT, [0,1], [0,1], XSTYLE=5, YSTYLE=5, COLOR=1, /NODATA
      XYOUTS, 0.1, 0.95, '!6n!DPOD!N = ' + rm0es(par.n_pod_triplets), $
        COLOR=1, /NORMAL
      XYOUTS, 0.1, 0.9, '!6n!Dfields!N = ' + rm0es(par.n_fields), $
        COLOR=1, /NORMAL

      n_triplets_thispage = p NE n_pages - 1 ? n_triplets_pp : $
        par.n_triplets - n_triplets_pp * (n_pages - 1)

      FOR j = 0, 4 DO BEGIN
        XYOUTS, xcoord[j], 0.85, header[j], /NORMAL, COLOR=1

        FOR k = 0, n_triplets_thispage - 1 DO $
          XYOUTS, xcoord[j], $
          0.8 * (1 - 1.0 * k / n_triplets_pp), $
          triplets_str[k,j], COLOR=1, /NORMAL
      ENDFOR
    ENDFOR
  ENDIF

  set_output, diag, /reset

END
