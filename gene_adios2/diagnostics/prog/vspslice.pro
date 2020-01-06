FUNCTION vspslice_info

  RETURN, {$
    type      : 'vsp',$
    title     : 'vsp slices',$
    help_text : ['Averages over vsp variable and '+$
                 'plots the result in 1D or 2D (or prints the '+$
                 'resulting 0D value)' ],$
    ext_vars  : [['var','0','variable index (0 = G_es, 1 = G_em, 2 = Q_es, 3 = Q_em, 4 = abs(f_)); '+$
                  'default: 4'],$
		 ['zind','0','z index; default: -2 (all); '+$
                  '-1 for average'],$
                 ['vpind','0','vp index; default: -2 (all); '+$
                  '-1 for average'],$
                 ['muind','0','mu index; default: nw0/2; -1 for '+$
                  'average; -2 for all']]}

END

;######################################################################

PRO vspslice_init, diag

  COMMON global_vars

  e = (*diag).external_vars

  IF N_ELEMENTS(*e.var) LT 1 THEN *e.var = 4
  IF N_ELEMENTS(*e.zind) LT 1 THEN *e.zind = -2
  IF N_ELEMENTS(*e.vpind) LT 1 THEN *e.vpind = -2
  IF N_ELEMENTS(*e.muind) LT 1 THEN *e.muind = par.nw0 / 2


  IF (*e.zind)[0] EQ -2 THEN *e.zind = INDGEN(par.nz0)
  IF (*e.vpind)[0] EQ -2 THEN *e.vpind = INDGEN(par.nv0)
  IF (*e.muind)[0] EQ -2 THEN *e.muind = INDGEN(par.nw0)

  i = set_internal_vars(diag,{$
    var      : *e.var,$
    zind     : *e.zind,$
    vpind     : *e.vpind,$
    muind     : *e.muind,$
    data_id  : PTR_NEW()})

END

;######################################################################

PRO vspslice_loop, diag

  COMMON global_vars

  i = (*diag).internal_vars

  zavg = (*i).zind[0] EQ -1
  vpavg = (*i).vpind[0] EQ -1
  muavg = (*i).muind[0] EQ -1

  zind = zavg ? INDGEN(par.nz0) : (*i).zind
  vpind = vpavg ? INDGEN(par.nv0) : (*i).vpind
  muind = muavg ? INDGEN(par.nw0) : (*i).muind

  zdim = N_ELEMENTS(zind)
  vpdim = N_ELEMENTS(vpind)
  mudim = N_ELEMENTS(muind)

  n_vars = N_ELEMENTS((*i).var)

  data = DBLARR(zavg?1:zdim,vpavg?1:vpdim,muavg?1:mudim,n_vars,gui.out.n_spec_sel)

  FOR isp = 0, gui.out.n_spec_sel - 1 DO BEGIN
    sp = (*gui.out.spec_select)[isp]

    FOR ivar = 0, n_vars - 1 DO BEGIN

       zdim = N_ELEMENTS(zind)
       vpdim = N_ELEMENTS(vpind)
       mudim = N_ELEMENTS(muind)
       ;datptr = vsp[isp,(*i).var[ivar]]

       ; select only required data
       ;data_sp = REFORM((*datptr)[zind,vpind,muind,0],$
       ;                 [zdim,vpdim,mudim],/OVERWRITE)
       ;data_sp = REFORM(time_avg((*i).data_id,(*vsp)[zind,vpind,muind,isp,$
       ;      (*i).var],vsp_time,fwd=gui.out.res_steps),$
       ;      [zdim,vpdim,mudim],/OVERWRITE)
       data_sp = REFORM((*vsp)[zind,vpind,muind,isp,$
             (*i).var],[zdim,vpdim,mudim],/OVERWRITE)

       ; now take absolute square or real/imaginary part
       ; note: Fourier averages are of quadratic type
       ;       regardless of the modifier choice
       ;data_sp = data_modifier(diag,data_sp)

       IF muavg THEN BEGIN
          ;todo more sphisticated integral!!
          GetMuWeightsAndKnots, mu_weights, mu
          data_sp = TEMPORARY(data_sp) * REBIN(REFORM($
                    mu_weights[muind,0,0],$
                    [1,1,mudim],/OVERWRITE),[zdim,vpdim,mudim])
          data_sp = REFORM(TOTAL(TEMPORARY(data_sp),3),$
               [zdim,vpdim,1],/OVERWRITE)
          mudim = 1
       ENDIF

       IF vpavg THEN BEGIN
          data_sp = REFORM(TOTAL(TEMPORARY(data_sp),2),$
                [zdim,1,mudim],/OVERWRITE) / par.nv0
          vpdim = 1
       ENDIF

       IF zavg THEN BEGIN
          data_sp = REFORM(TOTAL(TEMPORARY(data_sp),1)/par.nz0,$
                  [1,vpdim,mudim],/OVERWRITE)
          zdim = 1
       ENDIF

       data[*,*,*,ivar,isp] = data_sp[*,*,*]

    ENDFOR  ; --- ivar loop
  ENDFOR    ; --- species loop

  (*i).data_id = time_avg((*i).data_id,data,vsp_time,$
                              fwd=gui.out.res_steps)

END

;######################################################################

PRO vspslice_output, diag

  COMMON global_vars

  i = (*diag).internal_vars

  zavg = (*i).zind[0] EQ -1
  vpavg = (*i).vpind[0] EQ -1
  muavg = (*i).muind[0] EQ -1

  nzind = zavg ? 1 : N_ELEMENTS((*i).zind)
  nvpind = vpavg ? 1 : N_ELEMENTS((*i).vpind)
  nmuind = muavg ? 1 : N_ELEMENTS((*i).muind)

  n_res_steps = gui.out.res_steps * (series.step_count - 1) + 1

  ;IF (*i).ky2y THEN y_axis = -series.ly/2+INDGEN(par.ny0)*series.dy $
  ;  ELSE y_axis = (*series.ky)

  tmparr = [nzind,nvpind,nmuind]

  ndims = TOTAL(tmparr GT 1)

  vpar = - par.lv + INDGEN(par.nv0) * 2.0 * par.lv / (par.nv0 - 1.0)
  GetMuWeightsAndKnots, mu_weights, mu

  var_str=['G_es', 'G_em', 'Q_es', 'Q_em', 'abs(f_)']

  IF ndims GT 0 THEN BEGIN
    choice = WHERE(tmparr GT 1)

    my_label = STRARR(ndims)
    my_axis = PTRARR(ndims)

    FOR j = 0, ndims - 1 DO BEGIN
      CASE choice[j] OF
        0 : BEGIN
              my_axis[j] = PTR_NEW((*par.z)[(*i).zind])
              my_label[j] = 'z/L!Dref!N'
              ;print,  'z',*my_axis[j]
            END
        1 : BEGIN
              my_axis[j] = PTR_NEW(vpar[(*i).vpind])
              my_label[j] = '!6v!I!9#!6!N'
              ;print,  'vp',*my_axis[j]
            END
        2 : BEGIN
              my_axis[j] = PTR_NEW(mu[(*i).muind])
              my_label[j] = '!7l!6'
              ;print,  'mu', *my_axis[j]
            END
      ENDCASE
    ENDFOR
  ENDIF

  zval=*par.z
  info='z'
  IF ((*i).zind)[0] EQ -1 THEN info += ' av.' ELSE IF nzind EQ 1 THEN $
    info += ' ='+rm0es(zval[(*i).zind],prec=3)
  IF ((*i).vpind)[0] EQ -1 THEN info += ', !6v!I!9#!6!N av.' ELSE IF nvpind EQ 1 THEN $
    info += ', !6v!I!9#!6!N='+rm0es(vpar[(*i).vpind],prec=3)
  IF ((*i).muind)[0] EQ -1 THEN info += ', !7l!6 avg.' ELSE IF nmuind EQ 1 THEN $
    info += ', !7l!6='+rm0es(mu[(*i).muind],prec=3)

  IF ndims NE 0 THEN set_output, diag, /ps, coltable=33, charsize=csize
  nsp = gui.out.n_spec_sel
  nvars = N_ELEMENTS((*i).var)
  first = 1

  rdims = [nvars,nsp,n_res_steps]
  IF (nmuind GT 1) THEN rdims = [nmuind,[rdims]]
  IF (nvpind GT 1) THEN rdims = [nvpind,[rdims]]
  IF (nzind GT 1) THEN rdims = [nzind,[rdims]]

  data = REFORM(time_avg((*i).data_id,/avg,fwd=gui.out.res_steps,$
    tarr=time),rdims)

  FOR n = 0, n_res_steps - 1 DO BEGIN
   timestr = (gui.out.res_steps ? 'time = '+rm0es(time[n]) : 'time avg.')
   FOR isp = 0, nsp - 1 DO BEGIN
    sp = (*gui.out.spec_select)[isp]

    FOR ivar= 0 , nvars - 1 DO BEGIN

      ;mytitle = data_modifier(diag,get_var_string((*i).var[ivar],/fancy)+$
      ;           ((*i).var[ivar] LT par.n_fields ? '' : ',' + spec[sp].name),/legend)
      mytitle = var_str[(*i).var[ivar]]+' '+spec[sp].name

      CASE ndims OF
         0 : BEGIN
            IF (isp+ivar) EQ 0 THEN PRINT, timestr
            namestr = 'var'+rm0es(ivar)
            PRINT, STRING(namestr,format='(A-16)') + rm0es(data[ivar,isp,n])
         END
         1 : BEGIN
            IF (ivar EQ 0) THEN BEGIN
               collect_dat = REFORM(data[*,ivar,isp,n])
               collect_names = mytitle ;rm_idl_fmt(mytitle)
            ENDIF ELSE BEGIN
               collect_dat = [[collect_dat], [REFORM(data[*,ivar,isp,n])]]
               collect_names = [[collect_names],[mytitle]] ;[rm_idl_fmt(mytitle)]]
            ENDELSE

            IF (ivar EQ nvars-1) THEN BEGIN
               vspslice_plot, collect_dat, *my_axis[0],$
                  xtitle = my_label[0], info=info, title=collect_names,$
                  pos = [0.125,0.1+0.9/nsp*(nsp-1-isp),0.925,1.0-0.9/nsp*isp], $
                  no_erase=(isp GT 0)

               set_output, diag, commentlines=['species: '+spec[sp].name+', '+timestr+$
                    ', '+rm_idl_fmt(info)], header=[[rm_idl_fmt(my_label[0])],$
                   [rm_idl_fmt(collect_names)]], $
                   append=(n GT 0), dat=[[*my_axis[0]],[collect_dat]]
               plot_info_str, diag, time=(gui.out.res_steps ? time[t] : undef)
            ENDIF
	  END
      2 : BEGIN
            vspslice_plot, REFORM(data[*,*,ivar,isp,n]), *my_axis[0], $
              *my_axis[1], xtitle = my_label[0], ytitle = my_label[1],$
              pos = [0.125,0.1+0.9/(nvars*nsp)*(nvars*nsp-1-nvars*isp-ivar),0.925,1.0-0.9/(nvars*nsp)*(nvars*isp+ivar)], $
              no_erase=(isp+ivar GT 0), info=info, title=mytitle
            set_output, diag, commentlines=['species: '+spec[sp].name+', '+timestr,$
              '2D array '+rm_idl_fmt(mytitle)+'('+rm_idl_fmt(my_label[0])+$
              ','+rm_idl_fmt(my_label[1])+') at '+info],$
              append=(first EQ 0), dat=[[REFORM(data[*,*,ivar,isp,n])]]
            plot_info_str, diag, time=(gui.out.res_steps ? time[n] : undef)
	  END
      3 : BEGIN
            vspslice_plot, REFORM(data[*,*,*,ivar,isp,n]), *my_axis[0], $
              *my_axis[1], *my_axis[2], xtitle=my_label[0], ytitle=my_label[1],$
              info=info, title=mytitle
            set_output, diag, commentlines=['species: '+spec[sp].name+', '+timestr,$
              '3D array '+rm_idl_fmt(mytitle)+'('+rm_idl_fmt(my_label[0])+$
              ','+rm_idl_fmt(my_label[1])+','+rm_idl_fmt(my_label[2])+') at '+info],$
              append=(first EQ 0), dat=[[REFORM(data[*,*,*,ivar,isp,n])]]
            plot_info_str, diag, time=(gui.out.res_steps ? time[n] : undef)
	  END
      ELSE : PRINT, ndims, ' dimensions not supported'
     ENDCASE

     first = 0

     ENDFOR  ; --- ivar loop
   ENDFOR ; --- species loop
  ENDFOR

  IF ndims NE 0 THEN BEGIN
    PTR_FREE, my_axis

    set_output, diag, /reset
  ENDIF

END

;######################################################################

;same as in slice.pro
PRO vspslice_plot, dat, xaxis, yaxis, zaxis, xtitle=xtitle, ytitle=ytitle,$
  title=title, info=info, pos=pos, no_erase=no_erase, oplot=oplot, $
  norm = norm


  IF NOT KEYWORD_SET(xtitle) THEN xtitle = ''
  IF NOT KEYWORD_SET(ytitle) THEN ytitle = ''
  IF NOT KEYWORD_SET(title) THEN title = ''
  IF NOT KEYWORD_SET(info) THEN info = ''
  IF NOT KEYWORD_SET(pos) THEN pos = [0.0,0.0,1.0,1.0]
  IF N_ELEMENTS(no_erase) NE 1 THEN no_erase = 1

  COMMON global_vars

  datmax = MAX(dat,min=datmin)
  dpos_x = pos[2] - pos[0]
  dpos_y = pos[3] - pos[1]
  dy = 0.15

  color = INDGEN(128)+1
  color[0] = 4
  color[5] = 1

  IF N_ELEMENTS(dat) EQ 1 THEN BEGIN
    ; do nothing
  ENDIF ELSE IF N_ELEMENTS(yaxis) LT 1 THEN BEGIN
    LOADCT, 41, FILE='internal/colortable.tbl'

    ndims = SIZE(dat,/n_dimensions)
    IF ndims GT 1 THEN nvars = N_ELEMENTS(dat[0,*]) $
    ELSE nvars = 1

    mydat = FLTARR(N_ELEMENTS(xaxis),ndims)
    mydat[*,*] = dat ;copy data to allow for interference free modifications

    IF (ndims GT 1) THEN BEGIN
       mytitle = 'norm.: '
       FOR ivar = 0, nvars-1 DO BEGIN
          mydat[*,ivar] /= MAX(ABS(mydat[*,ivar]))
          mytitle += title[ivar]+' '
       ENDFOR
    ENDIF ELSE BEGIN
       nvars = 1
       mytitle = title
    ENDELSE

    datmax = MAX(mydat,min=datmin)

    PLOT, xaxis, mydat[*,0],$
          XSTYLE=1,YSTYLE=1, XTITLE=xtitle,YTITLE=ytitle,$
          YRANGE=[datmin,datmax],$
          TITLE=mytitle, color=1, /NODATA, NOERASE=no_erase,$
          POSITION=[pos[0]+0.01*dpos_x,pos[1]+dy*dpos_y,$
                    pos[0]+0.95*dpos_x,pos[3]-dy*dpos_y]

    FOR ivar = 0, nvars - 1 DO $
       OPLOT, xaxis, mydat[*,ivar], color=color[ivar]

    IF PRODUCT(!Y.CRANGE) LT 0 THEN OPLOT, !X.CRANGE, [0,0],color=1
  ENDIF ELSE IF N_ELEMENTS(zaxis) LT 1 THEN BEGIN
    LOADCT, 33, FILE='internal/colortable.tbl'
    lev_col = contour_levels([datmin,datmax],c_levels)
    c_levels = 20
    CONTOUR,dat[*,*], xaxis,yaxis, NOERASE=no_erase,$
      LEVELS=lev_col[*,0], C_COLORS=lev_col[*,1], /FILL,$
      XSTYLE=1,YSTYLE=1, XTITLE=xtitle,YTITLE=ytitle,$
      POSITION=[pos[0]+0.01*dpos_x,pos[1]+dy*dpos_y,$
        pos[0]+0.9*dpos_x,pos[3]-dy*dpos_y],$
      title=title
    plot_colorbar, lev_col, charsize=0.75*!P.CHARSIZE,prec=2,$
      POSITION=[pos[0]+0.95*dpos_x,pos[1]+dy*dpos_y,$
        pos[2],pos[3]-dy*dpos_y]

 ;  LOADCT, 41, FILE='internal/colortable.tbl'

 ;  SURFACE, dat[*,*],xaxis, yaxis, $
 ;    COLOR=1, /XSTYLE,/YSTYLE, $ ;/LEGO, $
 ;    POSITION=[0.1,0.1,0.95,0.475],$
 ;    YTITLE=ytitle, XTITLE=xtitle ;, TITLE=title
  ENDIF ELSE BEGIN
    LOADCT, 33, FILE='internal/colortable.tbl'
    print, 'starting slicer ...'
    datptr = PTR_NEW(dat)
    SLICER3,datptr,/MODAL
    PTR_FREE, datptr
  ENDELSE

  XYOUTS, 0.005, 0.025, info, COLOR=1, CHARSIZE=1.0, /NORMAL

END
