FUNCTION epar_info

  RETURN, {$
    type      : 'mom',$
    title     : 'Parallel E field',$
    help_text : ['Computes E_par from available contributions dPhi/dz, '+$
                 'B_perp-caused flutter from E_perp, and d/dt A_par term.'],$
    ext_vars  : [['n_levels','0','set to number of contour levels '+$
                  'for 2D contour plots (at outboard midplane); '+$
                  'default: 0 (plot of rms Epar as function of time)'],$
                 ['sparse_plot','0','plot contours only every n-th step '+$
                  '(use instead of usual sparse factor for higher '+$
                  'dA_par/dt accuracy); default: 1 (plot all steps)'],$
                 ['show_heating','1','show j_par E_par heating data '+$
                  '(for fully physical results, select all species)']]}

END

;######################################################################

PRO epar_init, diag

  COMMON global_vars

  IF NOT par.nonlinear OR NOT par.x_local THEN BEGIN
    PRINT, (*diag).name + ' error: requires local nonlinear runs'
    (*diag).selected = 0
    RETURN
  ENDIF

  show_static = par.nz0 GE 3
  with_em = par.n_fields GE 2
  IF NOT show_static AND NOT with_em THEN BEGIN
    PRINT, (*diag).name + ' skipped: neither es nor em contributions to E_par'
    (*diag).selected = 0
    RETURN
  ENDIF

  IF gui.out.n_spec_sel NE par.n_spec THEN PRINT, (*diag).name + $
    ' warning: not all species used for j_par evaluation!'

  e = (*diag).external_vars
  n_levels = N_ELEMENTS(*e.n_levels) EQ 1 ? *e.n_levels : 0
  sparse_plot = N_ELEMENTS(*e.sparse_plot) EQ 1 ? *e.sparse_plot : 1
  show_heating = KEYWORD_SET(*e.show_heating)

  i = set_internal_vars(diag,{$
    n_levels     : n_levels,$
    sparse_plot  : sparse_plot,$
    show_static  : show_static,$
    with_em      : with_em,$
    show_heating : show_heating,$
    time_id      : PTR_NEW(),$
    static_id    : PTR_NEW(),$
    flutter_id   : PTR_NEW(),$
    Apar_id      : PTR_NEW(),$
    jpar_id      : PTR_NEW()})

  IF show_static THEN fft_format, sxsy=0
  IF with_em THEN fft_format, kxky=[0,1], sxsy=1
  IF show_heating THEN fft_format, sxsy=par.n_fields+5

END

;######################################################################

PRO epar_loop, diag

  COMMON global_vars

  i = (*diag).internal_vars

  ; change value after ":" to plot contours at other z than nz0 / 2
  ; unless nz0 = 1, do not use boundary points (0 or par.nz0 - 1)
  zarr = (*i).n_levels EQ 0 ? INDGEN(par.nz0) : par.nz0 / 2
  nz = N_ELEMENTS(zarr)

  ; --- compute dPhi/dz contribution ---
  IF (*i).show_static THEN BEGIN
    IF (*i).n_levels GT 0 THEN BEGIN ; real space contours; nz = 1
      dPhi_dz = ((*mom[0,0].sxsy)[*,*,zarr+1] - $
        (*mom[0,0].sxsy)[*,*,zarr-1]) / REBIN($
        [(*par.z)[zarr+1]-(*par.z)[zarr-1]],[par.nx0,par.ny0])
    ENDIF ELSE BEGIN ; rms over all z points; nz = par.nz0
      dPhi_dz = pf_arr([par.nx0,par.ny0,nz])
      dPhi_dz[0,0,1] = ((*mom[0,0].sxsy)[*,*,2:*] - $
        (*mom[0,0].sxsy)[*,*,0:nz-3]) / REBIN(REFORM($
        (*par.z)[2:*]-(*par.z)[0:nz-3],[1,1,nz-2]),[par.nx0,par.ny0,nz-2])
      ; boundary points
      dPhi_dz[0,0,0] = $
        ((*mom[0,0].sxsy)[*,*,1] - (*mom[0,0].sxsy)[*,*,0]) / $
        REBIN([(*par.z)[1]-(*par.z)[0]],[par.nx0,par.ny0])
      dPhi_dz[0,0,nz-1] = $
        ((*mom[0,0].sxsy)[*,*,nz-1] - (*mom[0,0].sxsy)[*,*,nz-2]) / $
        REBIN([(*par.z)[nz-1]-(*par.z)[nz-2]],[par.nx0,par.ny0])
    ENDELSE

    static = - TEMPORARY(dPhi_dz)
    (*i).static_id = store_step((*i).static_id,static)
  ENDIF

  IF (*i).with_em THEN BEGIN ; electromagnetic contributions
    ; --- compute B_perp E_perp contribution ---
    ; divide A_par by B_0
    Apar_norm = (*mom[0,1].kxky)[*,*,zarr] * REBIN(REFORM($
      [1.0/(*series.geom).Bfield[zarr]],[1,1,nz]),[par.nkx0,par.nky0,nz])

    temp_sort = COMPLEXARR(par.nx0,par.ny0,nz,/NOZERO)
    IF NOT par.nky0 MOD 2 THEN temp_sort[*,par.nky0,*] = 0

    ; obtain real space B_x
    temp_kxky = COMPLEX(0,1) * REBIN(REFORM(*series.ky,[1,par.nky0]),$
      [par.nx0,par.nky0,nz]) * Apar_norm
    temp_sxky = FFT(TEMPORARY(temp_kxky),DIMENSION=1,/INVERSE,DOUBLE=0)
    FOR y = 0, par.nky0 - 1 DO temp_sort[0,y,0] = temp_sxky[*,y,*]
    FOR y = 1, par.nky0 - 1 DO $
      temp_sort[0,par.ny0-y,0] = CONJ(temp_sxky[*,y,*])
    B_x = FLOAT(FFT(temp_sort,DIMENSION=2,/INVERSE,DOUBLE=0))

    ; obtain real space B_y
    temp_kxky = - COMPLEX(0,1) * REBIN((*series.kx)[*,0],$
      [par.nx0,par.nky0,nz]) * Apar_norm
    temp_sxky = FFT(TEMPORARY(temp_kxky),DIMENSION=1,/INVERSE,DOUBLE=0)
    FOR y = 0, par.nky0 - 1 DO temp_sort[0,y,0] = temp_sxky[*,y,*]
    FOR y = 1, par.nky0 - 1 DO $
      temp_sort[0,par.ny0-y,0] = CONJ(temp_sxky[*,y,*])
    B_y = FLOAT(FFT(temp_sort,DIMENSION=2,/INVERSE,DOUBLE=0))

    ; obtain real space E_y
    temp_kxky = - COMPLEX(0,1) * REBIN(REFORM(*series.ky,[1,par.nky0]),$
      [par.nx0,par.nky0,nz]) * (*mom[0,0].kxky)[*,*,zarr]
    temp_sxky = FFT(TEMPORARY(temp_kxky),DIMENSION=1,/INVERSE,DOUBLE=0)
    FOR y = 0, par.nky0 - 1 DO temp_sort[0,y,0] = temp_sxky[*,y,*]
    FOR y = 1, par.nky0 - 1 DO $
      temp_sort[0,par.ny0-y,0] = CONJ(temp_sxky[*,y,*])
    E_y = FLOAT(FFT(temp_sort,DIMENSION=2,/INVERSE,DOUBLE=0))

    ; obtain real space E_x
    temp_kxky = - COMPLEX(0,1) * REBIN((*series.kx)[*,0],$
      [par.nx0,par.nky0,nz]) * (*mom[0,0].kxky)[*,*,zarr]
    temp_sxky = FFT(TEMPORARY(temp_kxky),DIMENSION=1,/INVERSE,DOUBLE=0)
    FOR y = 0, par.nky0 - 1 DO temp_sort[0,y,0] = temp_sxky[*,y,*]
    FOR y = 1, par.nky0 - 1 DO $
      temp_sort[0,par.ny0-y,0] = CONJ(temp_sxky[*,y,*])
    E_x = FLOAT(FFT(temp_sort,DIMENSION=2,/INVERSE,DOUBLE=0))

    flutter = TEMPORARY(B_x) * TEMPORARY(E_x) + $
      TEMPORARY(B_y) * TEMPORARY(E_y)
    (*i).flutter_id = store_step((*i).flutter_id,flutter)

    ; --- store A_par for dA_par/dt contribution ---
    A_par = (*mom[0,1].sxsy)[*,*,zarr]
    (*i).Apar_id = store_step((*i).Apar_id,A_par)
  ENDIF

  ; --- obtain, store j_par for heating ---
  IF (*i).show_heating THEN BEGIN
    j_par = pf_arr([par.nx0,par.ny0,nz],/zero)

    FOR n = 0, gui.out.n_spec_sel - 1 DO j_par += $
      spec[(*gui.out.spec_select)[n]].charge * $
      spec[(*gui.out.spec_select)[n]].dens * $
      (*mom[n,par.n_fields+5].sxsy)[*,*,zarr]

    (*i).jpar_id = store_step((*i).jpar_id,j_par)
  ENDIF

  (*i).time_id = store_step((*i).time_id,mom_time)

END

;######################################################################

PRO epar_output, diag

  COMMON global_vars

  i = (*diag).internal_vars

  n_steps = series.step_count
  nz = (*i).n_levels EQ 0 ? par.nz0 : 1
  time = store_step((*i).time_id,/get,/reg_array)

  E_par = pf_arr([par.nx0,par.ny0,nz,n_steps],/zero)

  IF (*i).show_static THEN BEGIN
    static = REFORM(store_step((*i).static_id,/get,/reg_array),$
      [par.nx0,par.ny0,nz,n_steps],/OVERWRITE)
    E_par += static
  ENDIF

  show_flutter = 0
  show_inductive = 0
  IF par.n_fields GE 2 THEN BEGIN
    show_flutter = 1
    flutter = REFORM(store_step((*i).flutter_id,/get,/reg_array),$
      [par.nx0,par.ny0,nz,n_steps],/OVERWRITE)
    E_par += flutter

    IF n_steps LT 2 THEN BEGIN
      PRINT, (*diag).name + $
        ' warning: only one time step, cannot evaluate dA_par/dt'
    ENDIF ELSE BEGIN
      show_inductive = 1

      ; --- evaluate dA_par/dt term ---
      A_par = REFORM(store_step((*i).Apar_id,/get,/reg_array),$
        [par.nx0,par.ny0,nz,n_steps],/OVERWRITE)

      dApar_dt = pf_arr([par.nx0,par.ny0,nz,n_steps])
      dApar_dt[0,0,0,1] = $
        (A_par[*,*,*,1:n_steps-1] - A_par[*,*,*,0:n_steps-2]) * $
        REBIN(REFORM(1.0/(time[1:n_steps-1]-time[0:n_steps-2]),$
        [1,1,1,n_steps-1]),[par.nx0,par.ny0,nz,n_steps-1])
      A_par = 0
      ; first time step: use same value as for second time step
      dApar_dt[0,0,0,0] = dApar_dt[*,*,*,1]

      inductive = - TEMPORARY(dApar_dt)
      E_par += inductive
    ENDELSE
  ENDIF

  IF (*i).show_heating THEN BEGIN
    j_par = store_step((*i).jpar_id,/get,/reg_array)
    j_spec = STRJOIN(spec[*gui.out.spec_select].name,',')

    heating = REFORM(TEMPORARY(j_par)*E_par,[par.nx0,par.ny0,nz,n_steps])
  ENDIF

  n_plots = ((*i).n_levels EQ 0 ? 1 : (1 + (*i).show_static + show_flutter + $
    show_inductive)) + (*i).show_heating
  IF (n_plots MOD 2 NE 0) AND (n_plots GT 1) THEN n_plots += 1
  n_plots_x = n_plots LE 3 ? 1 : 2
  n_plots_y = n_plots / n_plots_x

  set_output, diag, /ps, coltable=((*i).n_levels EQ 0 ? 41 : 47), $
    multi=[0,n_plots_x,n_plots_y]

  IF (*i).n_levels EQ 0 THEN BEGIN ; rms plots
    E_par = SQRT(TOTAL(TOTAL(TOTAL(TEMPORARY(E_par)^2,1),1)*REBIN($
      (*series.geom).jac_norm/(1L*par.nx0*par.ny0*nz),[nz,n_steps]),1))

    Epar_avg_id = PTR_NEW()
    FOR n = 0, n_steps - 1 DO $
      Epar_avg_id = time_avg(Epar_avg_id,E_par[n],time[n])
    Epar_avg = time_avg(Epar_avg_id,/avg)

    IF (*i).show_static THEN BEGIN
      static = SQRT(TOTAL(TOTAL(TOTAL(TEMPORARY(static)^2,1),1)*REBIN($
        (*series.geom).jac_norm/(1L*par.nx0*par.ny0*nz),[nz,n_steps]),1))

      static_avg_id = PTR_NEW()
      FOR n = 0, n_steps - 1 DO $
        static_avg_id = time_avg(static_avg_id,static[n],time[n])
      static_avg = time_avg(static_avg_id,/avg)
    ENDIF ELSE BEGIN
      static = 0 * E_par
      static_avg = 0
    ENDELSE

    IF show_flutter THEN BEGIN
      flutter = SQRT(TOTAL(TOTAL(TOTAL(TEMPORARY(flutter)^2,1),1)*REBIN($
        (*series.geom).jac_norm/(1L*par.nx0*par.ny0*nz),[nz,n_steps]),1))

      flutter_avg_id = PTR_NEW()
      FOR n = 0, n_steps - 1 DO $
        flutter_avg_id = time_avg(flutter_avg_id,flutter[n],time[n])
      flutter_avg = time_avg(flutter_avg_id,/avg)
    ENDIF ELSE BEGIN
      flutter = 0 * E_par
      flutter_avg = 0
    ENDELSE

    IF show_inductive THEN BEGIN
      inductive = SQRT(TOTAL(TOTAL(TOTAL(TEMPORARY(inductive)^2,1),1)*REBIN($
        (*series.geom).jac_norm/(1L*par.nx0*par.ny0*nz),[nz,n_steps]),1))

      inductive_avg_id = PTR_NEW()
      FOR n = 0, n_steps - 1 DO $
        inductive_avg_id = time_avg(inductive_avg_id,inductive[n],time[n])
      inductive_avg = time_avg(inductive_avg_id,/avg)
    ENDIF ELSE BEGIN
      inductive = 0 * E_par
      inductive_avg = 0
    ENDELSE

    PRINT, (*diag).name + ': <E_par,total^2>^(1/2) = ' + $
      rm0es(Epar_avg,prec=4)
    PRINT, (*diag).name + ': <E_par,static^2>^(1/2) = ' + $
      rm0es(static_avg,prec=4)
    PRINT, (*diag).name + ': <E_par,flutter^2>^(1/2) = ' + $
      rm0es(flutter_avg,prec=4)
    PRINT, (*diag).name + ': <E_par,inductive^2>^(1/2) = ' + $
      rm0es(inductive_avg,prec=4)

    yrange = [0,MAX([E_par,static,flutter,inductive])]

    PLOT, time, E_par, COLOR=1, /XSTYLE, YRANGE=yrange, /YSTYLE, TITLE='!6', $
      XTITLE='!6'+get_var_string(0,/time,/units), YTITLE='!6E!D!9#!6!N'
    OPLOT, !X.CRANGE, [1,1] * Epar_avg, COLOR=1, LINE=1
    OPLOT, time, static, COLOR=2
    OPLOT, !X.CRANGE, [1,1] * static_avg, COLOR=2, LINE=1
    OPLOT, time, flutter, COLOR=4
    OPLOT, !X.CRANGE, [1,1] * flutter_avg, COLOR=4, LINE=1
    OPLOT, time, inductive, COLOR=6
    OPLOT, !X.CRANGE, [1,1] * inductive_avg, COLOR=6, LINE=1

    !P.MULTI[0] = 1 + (*i).show_heating
    PLOT, [0,1], [0,1], /NODATA, XSTYLE=5, YSTYLE=5
    ypos = 1.05 - 0.5 * !D.Y_CH_SIZE / (!D.Y_VSIZE * (1 + (*i).show_heating))
    OPLOT, [0.05,0.1], [1,1] * 1.05, COLOR=1, /NOCLIP
    XYOUTS, 0.13, ypos, '!6total', COLOR=1, CHARSIZE=1
    OPLOT, [0.25,0.3], [1,1] * 1.05, COLOR=2, /NOCLIP
    XYOUTS, 0.32, ypos, '!6static', COLOR=2, CHARSIZE=1
    OPLOT, [0.45,0.5], [1,1] * 1.05, COLOR=4, /NOCLIP
    XYOUTS, 0.52, ypos, '!6flutter', COLOR=4, CHARSIZE=1
    OPLOT, [0.65,0.7], [1,1] * 1.05, COLOR=6, /NOCLIP
    XYOUTS, 0.72, ypos, '!6inductive', COLOR=6, CHARSIZE=1

    plot_info_str, diag

    IF (*i).show_heating THEN BEGIN
      heating = SQRT(TOTAL(TOTAL(TOTAL(TEMPORARY(heating)^2,1),1)*REBIN($
        (*series.geom).jac_norm/(1L*par.nx0*par.ny0*nz),[nz,n_steps]),1))

      heating_avg_id = PTR_NEW()
      FOR n = 0, n_steps - 1 DO $
        heating_avg_id = time_avg(heating_avg_id,heating[n],time[n])
      heating_avg = time_avg(heating_avg_id,/avg)
      PRINT, (*diag).name + ': <(j_par E_par)^2>^(1/2) = ' + $
        rm0es(heating_avg,prec=4)

      yrange = [0,MAX(heating)]
      PLOT, time, heating, COLOR=1, /XSTYLE, /YSTYLE, YRANGE=yrange, $
        TITLE='!6', XTITLE='!6'+get_var_string(0,/time,/units), $
        YTITLE='!6j!D!9#!6!NE!D!9#!6'
      OPLOT, !X.CRANGE, [1,1] * heating_avg, COLOR=1, LINE=1
    ENDIF
  ENDIF ELSE BEGIN ; contour plots
    xaxis = (FINDGEN(par.nx0) / par.nx0 - 0.5) * par.lx
    yaxis = (FINDGEN(par.ny0) / par.ny0 - 0.5) * par.ly

    FOR n = 0, n_steps - 1 DO BEGIN
      IF n MOD (*i).sparse_plot EQ 0 THEN BEGIN
        n_plots_remaining = n_plots

        lev_col = contour_levels(E_par[*,*,0,n],(*i).n_levels)
        CONTOUR, E_par[*,*,0,n], xaxis, yaxis, /FILL, XSTYLE=5, YSTYLE=5, $
          LEVELS=lev_col[*,0], C_COLORS=lev_col[*,1], /ISOTROPIC
        plot_colorbar, lev_col, orientation=1

        store_colors, store_rgb, new_ct=41
        !P.MULTI[0] = n_plots_remaining
        n_plots_remaining -= 1
        PLOT, xaxis, yaxis, /XSTYLE, /YSTYLE, COLOR=1, /NODATA, $
          /ISOTROPIC, XTITLE='!6x / !7q!6!Ds!N', YTITLE='!6y / !7q!6!Ds!N', $
          TITLE='!6E!D!9#!6,total!N'
        store_colors, store_rgb, /restore

        IF (*i).show_static THEN BEGIN
          lev_col = contour_levels(static[*,*,0,n],(*i).n_levels)
          CONTOUR, static[*,*,0,n], xaxis, yaxis, /FILL, XSTYLE=5, YSTYLE=5, $
            LEVELS=lev_col[*,0], C_COLORS=lev_col[*,1], /ISOTROPIC
          plot_colorbar, lev_col, orientation=1

          store_colors, store_rgb, new_ct=41
          !P.MULTI[0] = n_plots_remaining
          n_plots_remaining -= 1
          PLOT, xaxis, yaxis, /XSTYLE, /YSTYLE, COLOR=1, /NODATA, $
            /ISOTROPIC, XTITLE='!6x / !7q!6!Ds!N', YTITLE='!6y / !7q!6!Ds!N', $
            TITLE='!6E!D!9#!6,static!N'
          store_colors, store_rgb, /restore
        ENDIF

        IF show_flutter THEN BEGIN
          lev_col = contour_levels(flutter[*,*,0,n],(*i).n_levels)
          CONTOUR, flutter[*,*,0,n], xaxis, yaxis, /FILL, XSTYLE=5, YSTYLE=5, $
            LEVELS=lev_col[*,0], C_COLORS=lev_col[*,1], /ISOTROPIC
          plot_colorbar, lev_col, orientation=1

          store_colors, store_rgb, new_ct=41
          !P.MULTI[0] = n_plots_remaining
          n_plots_remaining -= 1
          PLOT, xaxis, yaxis, /XSTYLE, /YSTYLE, COLOR=1, /NODATA, $
            /ISOTROPIC, XTITLE='!6x / !7q!6!Ds!N', YTITLE='!6y / !7q!6!Ds!N', $
            TITLE='!6E!D!9#!6,flutter!N'
          store_colors, store_rgb, /restore
        ENDIF

        IF show_inductive THEN BEGIN
          lev_col = contour_levels(inductive[*,*,0,n],(*i).n_levels)
          CONTOUR, inductive[*,*,0,n], xaxis, yaxis, /FILL, XSTYLE=5, YSTYLE=5, $
            LEVELS=lev_col[*,0], C_COLORS=lev_col[*,1], /ISOTROPIC
          plot_colorbar, lev_col, orientation=1

          store_colors, store_rgb, new_ct=41
          !P.MULTI[0] = n_plots_remaining
          n_plots_remaining -= 1
          PLOT, xaxis, yaxis, /XSTYLE, /YSTYLE, COLOR=1, /NODATA, $
            /ISOTROPIC, XTITLE='!6x / !7q!6!Ds!N', YTITLE='!6y / !7q!6!Ds!N', $
            TITLE='!6E!D!9#!6,inductive!N'
          store_colors, store_rgb, /restore
        ENDIF

        IF (*i).show_heating THEN BEGIN
          lev_col = contour_levels(heating[*,*,0,n],(*i).n_levels)
          CONTOUR, heating[*,*,0,n], xaxis, yaxis, /FILL, XSTYLE=5, YSTYLE=5, $
            LEVELS=lev_col[*,0], C_COLORS=lev_col[*,1], /ISOTROPIC
          plot_colorbar, lev_col, orientation=1

          store_colors, store_rgb, new_ct=41
          !P.MULTI[0] = n_plots_remaining
          n_plots_remaining -= 1
          PLOT, xaxis, yaxis, /XSTYLE, /YSTYLE, COLOR=1, /NODATA, $
            /ISOTROPIC, XTITLE='!6x / !7q!6!Ds!N', YTITLE='!6y / !7q!6!Ds!N', $
            TITLE='!6j!D!9#!6'+j_spec+' !NE!D!9#!6!N'
          store_colors, store_rgb, /restore
        ENDIF

        !P.MULTI[0] = 0

        plot_info_str, diag, time=time[n]
      ENDIF
    ENDFOR
  ENDELSE

  set_output, diag, /reset

END
