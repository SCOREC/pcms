FUNCTION RZfilter_signal_info

  RETURN, {$
    type      : 'mom',$
    title     : 'RZ Gaussian filter',$
    help_text : ['Apply Gaussian spots to RZ data (for synth. diagnostics)'],$
    ext_vars  : [$
      ['var','0','variable to display; default: Tperp'],$
      ['rpos','0','array of radial spot positions'],$
      ['zpos','0','array of vertical spot positions'],$
      ['res','0','user defined resolutions: nz[, nx[, ny]]; default :'+$
                 'nz0>1280, nx0>200, ny0'],$
      ['tintrp','0','time interpolation factor; default: 5']]}

END

;######################################################################

PRO RZfilter_signal_init, diag

  COMMON global_vars

  e = (*diag).external_vars

  ; -------------------------------------------------------------------
  ; free parameters which are not included in the external vars
  torR = 1.0      ; major R normalization, default: 1.0 / GENE value
  phi_cut = 0.0   ; cut position, 0 to 2 pi, default: 0.0

  IF par.lx EQ 0 THEN BEGIN
    printerror, 'Skipping ' + (*diag).name + ': ky dependent lx'
    (*diag).selected = 0
    RETURN
  ENDIF

  IF (gui.out.res_steps EQ 0) THEN BEGIN
     printerror, 'Skipping ' + (*diag).name + ': time average not possible'
    (*diag).selected = 0
    RETURN
  ENDIF

  get_torus_pars, 2, torR, torZ, r0, $
                  round_q0, bigM, rho ;(see tools.pro)
  IF par.x_local AND (par.n0_global LT 0) THEN BEGIN
    IF ((bigM * round_q0) MOD 1 GT 1e-3) THEN BEGIN
      printerror, 'Skipping ' + (*diag).name + ': no proper bigM'
      (*diag).selected = 0
      RETURN
    ENDIF
  ENDIF

  ;apply correction factor if C_y NE r0/q0; this is currently only
  ;necessary for Miller geometry, and therefore only considered in the code
  ;below for local simulations
  Cyq0_r0=1.0
  if par.magn_geometry eq 'miller' then Cyq0_r0=(*series.geom).C_y*par.q0/r0


  ;default FWHMs of spots (may be overwritten by rpos, zpos choice)
  lr = 0.01
  lz = 0.01

  IF N_ELEMENTS(*e.rpos) EQ 1 THEN BEGIN
     CASE (*e.rpos) OF
        -1 : BEGIN
             (*e.rpos) = [2.10650,2.11150,2.12125,2.12625,$
                        2.13600,2.14100,2.15075,2.15575,2.16550,2.17050]
             lr = 0.01/4.
             END
        -2 : (*e.rpos) = [2.14138,2.15038,2.15938,$
                          2.16838,2.17738,2.14138,2.15038,2.15938,2.16838,2.17738,$
                          2.14138,2.15038,2.15938,2.16838,2.17738,2.14138,2.15038,$
                          2.15938,2.16838,2.17738,2.14138,2.15038,2.15938,2.16838,$
                          2.17738,2.14138,2.15038,2.15938,2.16838,2.17738]
        -3 : (*e.rpos) = [ 1.97444292,1.97053245,1.91387987,$
                           1.91166368,1.94373174,1.9399735 ]
        -4 : BEGIN
             (*e.rpos) = [2.04236079,2.03833697,1.98072508,$
                        1.97850208,2.01082238,2.00703401]
             lr = 0.0041627
             END
        -44 : BEGIN
             (*e.rpos) = [2.04236079,2.03833697,$
                        2.01082238,2.00703401]
             lr = 0.0041627
             END
        -5 : BEGIN
;            The radial spacing of the spot pair centres can be either 2mm or 4mm.
;            The spot pair separation and the array center is flexible. Here, we use
;            R(z=0) as center and a separation of 2.5cm (roughly as before)
             nspot = 6
             (*e.rpos) = (-(nspot-1.0)/2+FINDGEN(nspot))*0.025*0.5+(*series.geom).R[par.nz0/2]
             (*e.rpos)[1+2*INDGEN(nspot/2)] = (*e.rpos)[2*INDGEN(nspot/2)]+0.002
             ; Gaussian of 1/e radius of 2mm,  this is for 500eV and 2e19.
             ; The typical range is about 1.5-2.5mm.
             lr = 2*SQRT(ALOG(2.0)) * 0.002
             END
        -55 : BEGIN
;            The 1/e radial widths for this particular pair are 2.9mm.
;            The radial spacing is 4.2mm, this comes from the warm
;            resonance calculation, again including the filter bandwidths
;            Here, we use R(z=0) as center and a separation of 2.5cm (roughly as before)
             nspot = 6
             (*e.rpos) = (-(nspot-1.0)/2+FINDGEN(nspot))*0.025*0.5+(*series.geom).R[par.nz0/2]
             (*e.rpos)[1+2*INDGEN(nspot/2)] = (*e.rpos)[2*INDGEN(nspot/2)]+0.0042
             ; Gaussian of 1/e radius of 2mm,  this is for 500eV and 2e19.
             ; The typical range is about 1.5-2.5mm.
             lr = 2*SQRT(ALOG(2.0)) * 0.0029
          END
        -56 : BEGIN
;            Considering forward modeling as suggested by tbg,
;            including all the effect of the
;            filter bandwidth and a finite toroidal angle. The fitted
;            1/e width is now more like 4.7mm
;            The radial spacing is 4.2mm, this comes from the warm
;            resonance calculation, again including the filter bandwidths
;            Here, we use R(z=0) as center and a separation of 2.5cm (roughly as before)
             nspot = 6
             (*e.rpos) = (-(nspot-1.0)/2+FINDGEN(nspot))*0.025*0.5+(*series.geom).R[par.nz0/2]
             (*e.rpos)[1+2*INDGEN(nspot/2)] = (*e.rpos)[2*INDGEN(nspot/2)]+0.0042
             ; Gaussian of 1/e radius of 2mm,  this is for 500eV and 2e19.
             ; The typical range is about 1.5-2.5mm.
             lr = 2*SQRT(ALOG(2.0)) * 0.0047
          END
        -60 : BEGIN ;single pair
           nspot = 2
           (*e.rpos) = FLTARR(2)
           (*e.rpos)[0] = 2.061
           (*e.rpos)[1] = 2.065
           lr = 2*SQRT(ALOG(2.0)) * 0.0029
        END

        -6 : BEGIN
           (*e.rpos) = [ 2.12721474, 2.12476256, 2.12231585, 2.11988047, 2.11745039, 2.11502597,$
                         2.11233198, 2.10970354, 2.10708271, 2.10450591, 2.1019383, 2.09938135,$
                         2.09685995, 2.09434862, 2.09184396, 2.08936086, 2.08689032, 2.08437387,$
                         2.08217073, 2.07997963, 2.07779309, 2.07561149, 2.07344151, 2.07127847]
           lr = 2*SQRT(ALOG(2.0))*$
                [ 0.00152121, 0.00155259, 0.00147172, 0.00148436, 0.0015004, 0.00151622, 0.00152981,$
                  0.00154537, 0.00156143, 0.00157748, 0.00159037, 0.00160647,0.00162276, 0.00163784,$
                  0.00165216, 0.00166712, 0.00168799, 0.00170287, 0.00171948,0.00173791, 0.00175661,$
                  0.00177179, 0.0017911,  0.00181068]
        END
        -7 : BEGIN
           (*e.rpos) = [ 2.09084639, 2.08837207, 2.08590372, 2.08349036, 2.08129379, 2.07910442, $
                         2.07692014, 2.07474148, 2.07257591, 2.07041474, 2.06825855, 2.06611083, $
                         2.06397251, 2.06183908, 2.0597104, 2.05759227, 2.05548129, 2.05337534, $
                         2.05126838, 2.04917414, 2.0470849, 2.04500019, 2.04292167, 2.04085487]
           lr = 2*SQRT(ALOG(2.0))*$
                [ 0.00152121, 0.00155259, 0.00147172, 0.00148436, 0.0015004, 0.00151622, 0.00152981, $
                  0.00154537, 0.00156143, 0.00157748, 0.00159037, 0.00160647, 0.00162276, 0.00163784,$
                  0.00165216, 0.00166712, 0.00168799, 0.00170287, 0.00171948, 0.00173791, 0.00175661, $
                  0.00177179, 0.0017911, 0.00181068 ]
        END
        -8 : BEGIN
           (*e.rpos) = [ 2.05886125, 2.05674724, 2.05464208, 2.05253184, 2.05043, 2.04833798, 2.04625059,$
                         2.04416738, 2.04209436, 2.04002962, 2.0379688, 2.03591254, 2.0338675, 2.03182873,$
                         2.02979401, 2.02776378, 2.02574565, 2.02364947, 2.02151167, 2.01938106, $
                         2.01726135, 2.01514616, 2.01303547, 2.01093322]
           lr = 2*SQRT(ALOG(2.0))*$
                [ 0.00165848, 0.00167632, 0.00169464, 0.00170891, 0.00172695, 0.00174523, 0.00176293,$
                  0.00177936, 0.00179885, 0.00181876, 0.0018359, 0.0018545, 0.00187462, 0.00189473, $
                  0.00191072, 0.00196701, 0.00201861, 0.00203902, 0.00205605, 0.00207614, 0.00209576, $
                  0.00211231, 0.00212892, 0.00214673]
        END
        -9 : BEGIN
           nxx = par.nx0/2
           nxx += (nxx MOD 2)
           (*e.rpos) = (*series.geom).R[par.nz0/2]+(INDGEN(nxx)-nxx/2)*$
                       (*series.geom).c1[par.nz0/2]/(*series.geom).gxx[par.nz0/2]*par.lx/(nxx)*rho*par.Lref
           lr = 0.0001
        END
     ENDCASE
  ENDIF
  IF N_ELEMENTS(*e.zpos) EQ 1 THEN BEGIN
     CASE (*e.zpos) OF
        -1 : BEGIN
             (*e.zpos) = 0.049+FLTARR(10)
             (*e.zpos)[9] = 0.049001
             lz = 0.038/4. ;(in PoP17, (2010): 0.032/4 ?)
             END
        -2 : (*e.zpos) = [-0.04,-0.04,-0.04,-0.04,-0.04,$
                          -0.052,-0.052,-0.052,-0.052,-0.052,$
                          -0.064,-0.064,-0.064,-0.064,-0.064,$
                          -0.076,-0.076,-0.076,-0.076,-0.076,$
                          -0.088,-0.088,-0.088,-0.088,-0.088,$
                          -0.1,-0.1,-0.1,-0.1,-0.1]
        -3 : (*e.zpos) = [ 0.20057172,0.1998822,0.18989282,$
                           0.18950205,0.19515651,0.19449383]
        -4 : BEGIN
             (*e.zpos) = [0.21254747,0.21183796,0.20167943,$
                          0.20128746,0.2069864,0.20631841]
             lz = 0.0435642
             END
        -44 : BEGIN
             (*e.zpos) = [0.21254747,0.21183796,$
                          0.2069864,0.20631841]
             lz = 0.0435642
             END
        -5 : BEGIN
             ;poloidal steering angle; +4/-10 degrees most commonly used
;             theta = 4
             theta = -10 ;nT-phase must be done at -10.
             m = tan(-theta*!PI/180.)
             (*e.zpos) = m*(*e.rpos) + (0.301 - m*2.544)
             ;vertical 1/e^2 radius is about 4cm, it could be up to 20% bigger.
             ;With the new mirror we expect to have ~2cm at -10 degrees and ~1.6cm at +4 degrees.
             lz = SQRT(2*ALOG(2))*0.04
             END
        -55 : BEGIN
             ;poloidal steering angle; +4/-10 degrees most commonly used
;             theta = 4
             theta = -10 ;nT-phase must be done at -10.
             m = tan(-theta*!PI/180.)
             (*e.zpos) = m*(*e.rpos) + (0.301 - m*2.544)
             ;vertical 1/e^2 radius is about 4cm, it could be up to 20% bigger.
             ;With the new mirror we expect to have ~2cm at -10 degrees and ~1.6cm at +4 degrees.
             lz = SQRT(2*ALOG(2))*0.04*1.2
          END
        -56 : BEGIN
             ;poloidal steering angle; +4/-10 degrees most commonly used
;             theta = 4
             theta = -10 ;nT-phase must be done at -10.
             m = tan(-theta*!PI/180.)
             (*e.zpos) = m*(*e.rpos) + (0.301 - m*2.544)
             ;the latest revision of our beam spot size (1/e^2 radius) is 1.5 - 2.5 cm.
             lz = SQRT(2*ALOG(2))*0.015
           END
        -57 : BEGIN
             ;poloidal steering angle; +4/-10 degrees most commonly used
;             theta = 4
             theta = -10 ;nT-phase must be done at -10.
             m = tan(-theta*!PI/180.)
             (*e.zpos) = m*(*e.rpos) + (0.301 - m*2.544)
             ;the latest revision of our beam spot size (1/e^2 radius) is 1.5 - 2.5 cm.
             lz = SQRT(2*ALOG(2))*0.025
          END
        -58 : BEGIN
             ;poloidal steering angle; +4/-10 degrees most commonly used
;             theta = 4
             theta = -10 ;nT-phase must be done at -10.
             m = tan(-theta*!PI/180.)
             (*e.zpos) = m*(*e.rpos) + (0.301 - m*2.544)
             ;vertical 1/e^2 radius is about 4cm, it could be up to 20% bigger.
             ;here, we test 50% larger vertical radius
             lz = SQRT(2*ALOG(2))*0.04*1.5
          END

        -60 : BEGIN
           (*e.zpos) = FLTARR(2)
           (*e.zpos)[0] = 0.2158
           (*e.zpos)[1] = 0.2165
           lz = SQRT(2*ALOG(2))*0.015
        END
        -555 : BEGIN
             ;poloidal steering angle; +4/-10 degrees most commonly used
;             theta = 4
             theta = -10 ;nT-phase must be done at -10.
             m = tan(-theta*!PI/180.)
             (*e.zpos) = m*(*e.rpos) + (0.301 - m*2.544)
             ;vertical 1/e^2 radius is about 4cm, it could be up to 20% bigger.
             ;With the new mirror we expect to have ~2cm at -10 degrees and ~1.6cm at +4 degrees.
             lz = SQRT(2*ALOG(2))*0.02
          END
        -6 : BEGIN
           (*e.zpos) = [ 0.22750951, 0.22707713, 0.22664571, 0.22621628, 0.2257878, 0.2253603, $
                         0.22488528, 0.22442182, 0.22395969, 0.22350533, 0.22305259, 0.22260174, $
                         0.22215714, 0.22171433, 0.22127269, 0.22083485, 0.22039923, 0.21995551, $
                         0.21956704, 0.21918069, 0.21879514, 0.21841047, 0.21802784, 0.21764644]
           lz = SQRT(2*ALOG(2))*0.02
        END
        -7 : BEGIN
           (*e.zpos) = [ 0.22109679, 0.2206605, 0.22022527, 0.21979972, 0.21941241, 0.21902637, $
                         0.21864122, 0.21825706, 0.21787521, 0.21749414, 0.21711395, 0.21673525,$
                         0.2163582, 0.21598202, 0.21560668, 0.21523319, 0.21486097, 0.21448964, $
                         0.21411812, 0.21374885, 0.21338046, 0.21301287, 0.21264637, 0.21228194]
           lz = SQRT(2*ALOG(2))*0.02
        END
        -77 : BEGIN
           (*e.zpos) = [ 0.22109679, 0.2206605, 0.22022527, 0.21979972, 0.21941241, 0.21902637, $
                         0.21864122, 0.21825706, 0.21787521, 0.21749414, 0.21711395, 0.21673525,$
                         0.2163582, 0.21598202, 0.21560668, 0.21523319, 0.21486097, 0.21448964, $
                         0.21411812, 0.21374885, 0.21338046, 0.21301287, 0.21264637, 0.21228194]
           lz = SQRT(2*ALOG(2))*0.04
        END
        -8 : BEGIN
           (*e.zpos) = [ 0.21545695, 0.21508419, 0.214713, 0.2143409, 0.21397029, 0.21360141, $
                         0.21323335, 0.21286602, 0.21250049, 0.21213642, 0.21177305, 0.21141047, $
                         0.21104988, 0.21069039, 0.21033161, 0.20997363, 0.20961778, 0.20924816, $
                         0.20887121, 0.20849553, 0.20812176, 0.2077488, 0.20737663, 0.20700594]
           lz = SQRT(2*ALOG(2))*0.02
        END
        -88 : BEGIN
           (*e.zpos) = [ 0.21545695, 0.21508419, 0.214713, 0.2143409, 0.21397029, 0.21360141, $
                         0.21323335, 0.21286602, 0.21250049, 0.21213642, 0.21177305, 0.21141047, $
                         0.21104988, 0.21069039, 0.21033161, 0.20997363, 0.20961778, 0.20924816, $
                         0.20887121, 0.20849553, 0.20812176, 0.2077488, 0.20737663, 0.20700594]
           lz = SQRT(2*ALOG(2))*0.04
        END
         -9 : BEGIN
           (*e.zpos) = INTARR(N_ELEMENTS((*e.rpos)))+torZ
           lz = 0.0001 ;SQRT(2*ALOG(2))*0.02
        END
     ENDCASE
  ENDIF

  IF N_ELEMENTS(*e.rpos) NE N_ELEMENTS(*e.zpos) THEN BEGIN
    printerror, 'Skipping ' + (*diag).name + ': R/Z coordinate mismatch'
    (*diag).selected = 0
    RETURN
  ENDIF


  IF NOT KEYWORD_SET(*e.var) THEN *e.var = par.n_fields+2
  res = [par.nx0>200,par.ny0,par.nz0>1280]

  CASE N_ELEMENTS(*e.res) OF
    1    : res[2] = *e.res
    2    : res[*] = [(*e.res)[1],res[1],(*e.res)[0]]
    3    : res[*] = [(*e.res)[1],(*e.res)[2],(*e.res)[0]]
    ELSE :
  ENDCASE

  IF (N_ELEMENTS(*e.tintrp) LT 1) THEN *e.tintrp = 5


  use_man_minmax = N_ELEMENTS(man_minmax) EQ 2

  i = set_internal_vars(diag,{$
    var            : *e.var,$
    tintrp         : *e.tintrp,$
    x_res          : res[0],$
    y_res          : res[1],$
    z_res          : res[2],$
    n_vars         : N_ELEMENTS(*e.var),$
    torR           : torR,$
    torZ           : torZ,$
    r0      	   : r0,$
    rpos           : *e.rpos,$
    zpos           : *e.zpos,$
    lr             : lr,$
    lz             : lz,$
    dx             : 0.0,$
    dy             : 0.0,$
    q0             : round_q0,$
    bigM           : bigM,$
    y_shift        : get_y_shift(res[0],0.01, 0.01, Cyq0_r0),$
    closest_y_cut  : FLTARR(2,res[2],res[0],/NOZERO),$
    dy_cut         : FLTARR(2,res[2],res[0],/NOZERO),$
    coords         : get_RZ_coords(res[2], res[0], INDGEN(res[0]),$
                       torR=torR, torZ=torZ),$
    dat3d_id       : PTR_NEW(),$
    rho            : rho})

  use_sxky = 1

  (*i).dx = rho * par.lx / (*i).x_res
  (*i).dy = rho * par.ly / (*i).y_res
  (*i).y_shift = get_y_shift((use_sxky?par.nx0:(*i).x_res), $
       (*i).dx, (*i).dy, Cyq0_r0)

  ; calculating cut coordinates
  x_axis = (INDGEN((*i).x_res) - (*i).x_res / 2) * (*i).dx
  (*i).coords = get_RZ_coords((*i).z_res, (*i).x_res, x_axis,$
                  torR=(*i).torR, torZ=(*i).torZ)

  ; calculating y interpolation for cuts
  two_o_z_res = 2.0 / FLOAT((*i).z_res)
  dy_inv = 1.0 / (*i).dy
  yres_o_2 = (*i).y_res / 2

  FOR x = 0, (*i).x_res - 1 DO BEGIN
     IF par.x_local THEN BEGIN
        phi_auxterm = phi_cut * (*series.geom).C_y
        qprof_auxterm = !PI * ((*i).r0 + par.shat * $
                               x_axis[x])*Cyq0_r0
     ENDIF ELSE BEGIN
        phi_auxterm = phi_cut * C_y[x]
        qprof_auxterm = !PI * q_prof[x] * C_y[x]
     ENDELSE
     qprof_auxterm *= par.sign_Ip_CW

     slice_y = (qprof_auxterm * (two_o_z_res * INDGEN((*i).z_res) - $
                                 1.0) - phi_auxterm) * dy_inv + yres_o_2
     slice_y = ((slice_y MOD (*i).y_res) + (*i).y_res) MOD (*i).y_res

     closest_y = FIX(slice_y)
     dy_lo = slice_y - closest_y
     dy_hi = 1.0 - dy_lo
     closest_y_hi = (((closest_y + 1) MOD (*i).y_res) + $
                     (*i).y_res) MOD (*i).y_res
     closest_y_lo = ((closest_y MOD (*i).y_res) + $
                     (*i).y_res) MOD (*i).y_res

     (*i).closest_y_cut[0,*,x] = closest_y_lo
     (*i).closest_y_cut[1,*,x] = closest_y_hi
     (*i).dy_cut[0,*,x] = dy_hi
     (*i).dy_cut[1,*,x] = dy_lo
  ENDFOR

  fft_format, sxky=use_sxky*(*e.var), sxsy=(use_sxky EQ 0)*(*e.var)

END

;######################################################################

PRO RZfilter_signal_loop, diag
; first loop just storing data
; NO DOPPLER SHIFT SHOULD BE APPLIED SINCE THIS WILL BE DONE LATER
; WITH THE TIME INTERPOLATED DATA!!

  COMMON global_vars

  i = (*diag).internal_vars

  dat3d = COMPLEXARR(par.nx0,par.nky0,par.nz0,gui.out.n_spec_sel,/NOZERO)

  FOR isp = 0, gui.out.n_spec_sel - 1 DO BEGIN
     dat3d[*,*,*,isp] = (*mom[isp,(*i).var].sxky)[*,*,*]
     ;add white noise
;     noise = RANDOMN(SEED,par.nx0,par.nky0,par.nz0)
;     noise /= MAX(ABS(noise))
;     dat3d[*,*,*,isp] *= (1.0+0.5*exp(complex(0.0,2.0*!PI*noise)))
  ENDFOR

  (*i).dat3d_id = time_avg((*i).dat3d_id,dat3d,mom_time,fwd=gui.out.res_steps)

END
;######################################################################

PRO RZfilter_signal_output, diag

  COMMON global_vars

  IF (series.step_count LT 2) THEN RETURN

  i = (*diag).internal_vars

  NR = N_ELEMENTS((*i).rpos)
  NZ = N_ELEMENTS((*i).zpos)
  nsigs = NR

  ;get and time-interpolate fluxtube data
  datsxky_tintrp = time_avg((*i).dat3d_id,/avg,fwd=gui.out.res_steps,$
      tarr=time,intrp=(*i).tintrp)

  ntime = N_ELEMENTS(time)
  datsxky_tintrp = REFORM(datsxky_tintrp,[par.nx0,par.nky0,par.nz0,$
                  gui.out.n_spec_sel,ntime],/OVERWRITE)

  ;CHECK: do we need to apply the phase factor to kxky??
  IF (par.omegatorref NE 0.0) THEN BEGIN
     print, 'Applying background rotation ...'
;    NOTE: changed sign in June 2018
     phasefac = DCOMPLEX(0.0,1.0)*(*series.geom).C_y*$
             par.Lref/par.rhoref*par.Lref/par.cref
     IF (par.sign_Omega_CW NE 0) THEN $
        phasefac *= par.sign_Omega_CW*ABS(par.omegatorref) $
     ELSE phasefac *= par.omegatorref
     ;Attention: may give different results compared to standard ExB
     ;shift application on mom data since ExB_stime is not adapted
     ;in follow-up runs!
;ATTENTION: decided to ignore ExB_stime and always apply background
;rotation - otherwise, analysis of gamma_ExB=0 sims makes no sense
     FOR t = 0, ntime-1 DO BEGIN
;        IF time[t] GT par.ExB_stime THEN BEGIN
           FOR iky = 0, par.nky0 - 1 DO datsxky_tintrp[*,iky,*,*,t] *= $
;              EXP(phasefac*(*par.ky)[iky]*(time[t]-par.ExB_stime))
              EXP(phasefac*(*par.ky)[iky]*(time[t])) ;-par.ExB_stime))
;        ENDIF
     ENDFOR
  ENDIF


  ;-----------------------------------------------------------------------------------

  rsc = (par.rhoref / par.Lref)

  time_si = (time-time[0])*par.Lref/par.cref
  dt_min = min(time_si[1:*]-time_si[0:ntime-2],max=dt_max)
  dt_min_avg = TOTAL(time_si[1:*]-time_si[0:ntime-2])/(ntime-1.0)
  print, 'orig. data: min/max. freq. [kHz] ',rm0es(1.0E-3/time_si[ntime-1])+$
         '/'+rm0es(1.0E-3/dt_min)
  freq_nyquist = 0.5E-3/dt_max
  print, 'orig. data: Nyquist freq. [kHz] = '+rm0es(freq_nyquist)+'-'+rm0es(0.5E-3/dt_min)

  ;kept these to use the old code from CECE_signal.pro
  ntintrp = 1
  nteqi = ntime
  dt_eqi = time_si[ntime-1]/(nteqi-1.0)
  time_eqi = time_si

  nposf = nteqi/2
  dfreq = 1.0/time_eqi[nteqi-1]

  NWIND = FIX(2E3/dfreq)>5 ;FIX(22E3/dfreq)>1
  NFFT = nteqi/NWIND ;128<nteqi  ;nteqi/NWIND
  NWIND = nteqi/NFFT
  print, 'number of Hanning windows:  ', rm0es(NWIND)
  print, 'number of freq. in windows: ', rm0es(NFFT)
  dffreq = dfreq*NWIND
  print, 'df in Hanning windows and after interpolation [kHz]: ', rm0es(dffreq*1E-3)

  ;-----------------------------------------------------------------------------------

  print, ''
  print, 'Spatial interpolation, FFT and RZ filter ...'

  cut_data = FLTARR((*i).z_res,(*i).x_res,/NOZERO)
  filter_signal = FLTARR(nsigs,ntime,/NOZERO)
  raw_signal = FLTARR(nsigs,ntime,/NOZERO)
  raw_data = FLTARR(nsigs,/NOZERO)

  ; for storing cut data ...
  xind = REBIN(REFORM(INDGEN((*i).x_res),[1,(*i).x_res]),$
               [(*i).z_res,(*i).x_res])
  zind = REBIN(INDGEN((*i).z_res),[(*i).z_res,(*i).x_res])
  cyc0 = REFORM((*i).closest_y_cut[0,*,*])
  cyc1 = REFORM((*i).closest_y_cut[1,*,*])

  FOR isp = 0, gui.out.n_spec_sel - 1 DO BEGIN
     sp=(*gui.out.spec_select)[isp]
     prev = 0

     set_output, diag,  sp, coltable=33

     FOR t = 0, ntime-1 DO BEGIN
        data = REFORM(datsxky_tintrp[*,*,*,isp,t])
;print, 'interpol ...'
        interp_3d_data, data, (*i).x_res, (*i).y_res, (*i).z_res,$
                  (par.x_local?(*i).y_shift/(*i).y_res:$
                   (*i).bigM*(*series.geom).q),infft=[0,1,0]
;print, 'computing cut_data ...'
        cut_data[0,0] = REFORM((*i).dy_cut[0,*,*]*$
           data[xind,cyc0,zind]+(*i).dy_cut[1,*,*]*data[xind,cyc1,zind])

        IF (t EQ 0) THEN BEGIN
           lev_col = contour_levels(cut_data,n_levels)
           CONTOUR, cut_data, (*i).coords.R,(*i).coords.Z, $
                    C_COLORS=lev_col[*,1], $
                    LEVELS=lev_col[*,0], /FILL, /ISOTROPIC, $
                    XSTYLE=1, YSTYLE=1
           plot_info_str, diag, time=time[t]
        ENDIF
;print, 'applying filter ...'
        RZ_arb_filter, cut_data, (*i).coords.R, (*i).coords.Z, $
            rpos=(*i).rpos,zpos=(*i).zpos,lr=(*i).lr,lz=(*i).lz,$
            raw_sig=raw_data, sig_data = filter_sig, oplot=(t EQ 0),$
            plot=(t EQ 0)
;print, 'assigning signals ...'
        filter_signal[*,t] = filter_sig ;synthetic data considering PSFs
        raw_signal[*,t] = raw_data  ;raw GENE data at spot locations

        IF (t GT prev*ntime/10) THEN BEGIN
           print, rm0es(10*prev)+'%'
           prev += 1
        ENDIF
     ENDFOR

     raw_signal *= rsc
     filter_signal *= rsc

     ;-------------------------------------------------------------------
     freq_spect = 1
     IF (freq_spect) THEN BEGIN

        LOADCT, 41, FILE='internal/colortable.tbl'

        raw_aspect_avg = FLTARR(nfft/2+1)
        filter_aspect_avg = FLTARR(nfft/2+1)

        raw_xspect_avg = FLTARR(nfft/2+1)
        filter_xspect_avg = FLTARR(nfft/2+1)

        FOR isig = 0,nsigs-1 DO BEGIN
           raw_aspect_avg += CALC_RFSPECT1D(raw_signal[isig,*],$
                                            FREQ=ffreq, NFFT=NFFT)
           filter_aspect_avg += CALC_RFSPECT1D(filter_signal[isig,*],$
                                               FREQ=ffreq, NFFT=NFFT)
        ENDFOR

                                ;assuming double-spots!!!
        FOR ir = 0, nsigs/2-1 DO BEGIN
           raw_xspect_avg += CALC_RFSPECT1D(raw_signal[2*ir,*], $
                                            raw_signal[2*ir+1,*], $
                                            FREQ=ffreq, NFFT=NFFT) ;, ERR=err, $
;                                     COHERENCY=coh)
           filter_xspect_avg += CALC_RFSPECT1D(filter_signal[2*ir,*], $
                                               filter_signal[2*ir+1,*], $
                                               FREQ=ffreq, NFFT=NFFT) ;, ERR=err, $
;                                     COHERENCY=coh)
        ENDFOR

        dffreq *= 1E-3          ;translate from Hz to kHz
        time_eqi *= 1E3         ;translate from s to ms
        ffreq *= dffreq

        fidx = WHERE((ffreq GE 40) AND ffreq LE 400) ;freq. window of expt. diag.

        filter_sig_avg = TOTAL(filter_signal[*,*],1)/nsigs

        raw_aspect_avg /= nsigs*dffreq
        filter_aspect_avg /= nsigs*dffreq

        raw_xspect_avg /= (nsigs/2)*dffreq
        filter_xspect_avg /= (nsigs/2)*dffreq

        ; time traces
        ymin=MIN(filter_signal[*,*],max=ymax)
        PLOT, time_eqi, REFORM(filter_signal[0,*]),COLOR=1,/NODATA,$
              /XSTYLE,/YSTYLE, YRANGE=[ymin<0,ymax], $
              xtitle='t / ms',ytitle='fluc'
        FOR isig = 1, nsigs-1 DO $
           OPLOT, time_eqi, filter_signal[isig,*], COLOR=1+isig

        ; auto + cross spectra
        yrsc = 1E6
        ymin=MIN([raw_aspect_avg,filter_aspect_avg]*yrsc,max=ymax)

        PLOT, ffreq, raw_aspect_avg*yrsc, COLOR=1,/NODATA,$
              /XSTYLE,/YSTYLE, YRANGE=[ymin<0,ymax], $
              xtitle='f / kHz',ytitle='pow. spectra x 10!E-6!N / kHz',$
              POSITION=[0.15,0.231,0.95,0.95]
        OPLOT, ffreq, raw_aspect_avg*yrsc, COLOR=2, PSYM=-4
        OPLOT, ffreq, raw_xspect_avg*yrsc, COLOR=3, PSYM=-5
        OPLOT, ffreq, filter_xspect_avg*yrsc, COLOR=4, PSYM=-6
        plot_legend, INDGEN(3)+2, ['raw autopow.', 'raw crosspow.', 'synth. crosspow'], $
                     x_offset=[0.15,0.0,0.95,0.12], per_line = 2, PSYM=-4-INDGEN(3)

        PLOT, ffreq, raw_aspect_avg*yrsc, COLOR=1,/NODATA,$
              /XSTYLE,/YSTYLE, YRANGE=[ymin<0,ymax], XRANGE=[0.0,500.0],$
              xtitle='f / kHz',ytitle='pow. spectra x 10!E-6!N / kHz',$
              POSITION=[0.15,0.231,0.95,0.95]
        OPLOT, ffreq, raw_aspect_avg*yrsc, COLOR=2, PSYM=-4
        OPLOT, ffreq, raw_xspect_avg*yrsc, COLOR=3, PSYM=-5
        OPLOT, ffreq, filter_xspect_avg*yrsc, COLOR=4, PSYM=-6
        plot_legend, INDGEN(3)+2, ['raw autopow.', 'raw crosspow.', 'synth. crosspow'], $
                     x_offset=[0.15,0.0,0.95,0.12], per_line = 2, PSYM=-4-INDGEN(3)
     ENDIF

     ;radial correlation function and length
     radcorr = 0
     IF (radcorr) THEN BEGIN
        lag = INDGEN(2*nsigs-3)-(nsigs-2)
        t = 0
        rawcorr = A_CORRELATE(raw_signal[*,t],lag)
        filtercorr = A_CORRELATE(filter_signal[*,t],lag)
        FOR t = 1, ntime-1 DO BEGIN
           rawcorr += A_CORRELATE(raw_signal[*,t],lag)
           filtercorr += A_CORRELATE(filter_signal[*,t],lag)
        ENDFOR
        rawcorr /= ntime
        filtercorr /= ntime

        ;meaningful x - axis
        dr_avg = MEAN(ABS((*i).rpos[1:nsigs-1]-(*i).rpos[0:nsigs-2]))*1E3
        lagr = lag*dr_avg       ;/(2.0*nsigs-3.0)

        ;correlation length at half max
        amp = 0.5

        rawr1 = INTERPOL(lagr[0:nsigs-2],rawcorr[0:nsigs-2],amp)
        rawr2 = INTERPOL(lagr[nsigs-2:2*nsigs-4],rawcorr[nsigs-2:2*nsigs-4],amp)
        raw_corrfwhm = rawr2-rawr1
        print, 'raw correlation FWHM [mm]: ',raw_corrfwhm

        filterr1 = INTERPOL(lagr[0:nsigs-2],filtercorr[0:nsigs-2],amp)
        filterr2 = INTERPOL(lagr[nsigs-2:2*nsigs-4],filtercorr[nsigs-2:2*nsigs-4],amp)
        synth_corrfwhm = filterr2-filterr1
        print, 'synth correlation FWHM [mm]: ',synth_corrfwhm

        amp = 1.0/EXP(1.0)
        rawr1 = INTERPOL(lagr[0:nsigs-2],rawcorr[0:nsigs-2],amp)
        rawr2 = INTERPOL(lagr[nsigs-2:2*nsigs-4],rawcorr[nsigs-2:2*nsigs-4],amp)
        raw_corrinve = rawr2-rawr1
        print, 'raw correlation 1/e [mm]: ',raw_corrinve/2

        filterr1 = INTERPOL(lagr[0:nsigs-2],filtercorr[0:nsigs-2],amp)
        filterr2 = INTERPOL(lagr[nsigs-2:2*nsigs-4],filtercorr[nsigs-2:2*nsigs-4],amp)
        synth_corrinve = filterr2-filterr1
        print, 'synth correlation 1/e [mm]: ',synth_corrinve/2

        ;PLOT correlation functions
        ymin = MIN([rawcorr,filtercorr],max=ymax)
        PLOT, lagr, rawcorr, COLOR=1,$
              /XSTYLE, /YSTYLE, xtitle='!7D!6r!Iavg!N / mm',$
              ytitle = 'auto correlation', YRANGE=[ymin,ymax],$
              POSITION=[0.15,0.231,0.95,0.95]
        OPLOT, lagr, filtercorr, COLOR=2

        OPLOT, [rawr1,rawr1],[!Y.CRANGE[0],amp],COLOR=1,LINESTYLE=1
        OPLOT, [rawr2,rawr2],[!Y.CRANGE[0],amp],COLOR=1,LINESTYLE=1
        OPLOT, [-raw_corrinve/2,raw_corrinve/2],[amp,amp],COLOR=1,LINESTYLE=1
        OPLOT, [filterr1,filterr1],[!Y.CRANGE[0],amp],COLOR=2,LINESTYLE=2
        OPLOT, [filterr2,filterr2],[!Y.CRANGE[0],amp],COLOR=2,LINESTYLE=2
        OPLOT, [-synth_corrinve/2,synth_corrinve/2],[amp,amp],COLOR=2,LINESTYLE=2

        plot_legend, INDGEN(2)+1, ['raw corr.', 'synth. corr'], $
          x_offset=[0.15,0.0,0.95,0.12], per_line = 2 ;, PSYM=-4-INDGEN(3)

     ENDIF

     resstr = 'Resolution (nx,ny,nz,nt) = '+rm0es((*i).x_res)+', '+$
              rm0es((*i).y_res)+', '+rm0es((*i).z_res)+', '+rm0es((*i).tintrp)

     IF N_ELEMENTS((*i).rpos) GT 1 THEN BEGIN
        rposstr = 'rpos = '+rm0es((*i).rpos[0])
        FOR isig = 1, nsigs-1 DO rposstr += ', '+rm0es((*i).rpos[isig])
     ENDIF ELSE rposstr = 'rpos = '+rm0es((*i).rpos)
     IF N_ELEMENTS((*i).zpos) GT 1 THEN BEGIN
        zposstr = 'zpos = '+rm0es((*i).zpos[0])
        FOR isig = 1, nsigs-1 DO zposstr += ', '+rm0es((*i).zpos[isig])
     ENDIF ELSE zposstr = 'zpos = '+rm0es((*i).zpos)

     IF N_ELEMENTS((*i).lr) GT 1 THEN BEGIN
        lrstr = 'Lr = '+rm0es((*i).lr[0])
        FOR isig = 1, nsigs-1 DO lrstr += ', '+rm0es((*i).lr[isig])
     ENDIF ELSE lrstr = 'Lr = '+rm0es((*i).lr)
     IF N_ELEMENTS((*i).lz) GT 1 THEN BEGIN
        lzstr = 'Lz = '+rm0es((*i).lz[0])
        FOR isig = 1, nsigs-1 DO lzstr += ', '+rm0es((*i).lz[isig])
     ENDIF ELSE lzstr = 'Lz = '+rm0es((*i).lz)

     IF (freq_spect) THEN BEGIN
     infostr = [resstr, rposstr, zposstr, lrstr, lzstr,$
                'RZfilter - '+spec[sp].name+': sqrt(<|'+get_var_string((*i).var)+$
                '|^2>)/'+get_var_string((*i).var)+'_0:',$
               'raw/auto power/all freq.:    '+$
               rm0es(SQRT(INT_TABULATED(ffreq,raw_aspect_avg))*100,prec=3)+' % ',$
               'raw/cross power/all freq:    '+$
               rm0es(SQRT(INT_TABULATED(ffreq,raw_xspect_avg))*100,prec=3)+' %',$
               'raw/cross power/40-400kHz:   '+((fidx[0] LE 1)?'NaN':$
               rm0es(SQRT(INT_TABULATED(ffreq[fidx],raw_xspect_avg[fidx]))*100,prec=3)+' %'),$
               'synth/auto power/all freq.:  '+$
               rm0es(SQRT(INT_TABULATED(ffreq,filter_aspect_avg))*100,prec=3)+' %',$
               'synth/cross power/all freq.: '+$
               rm0es(SQRT(INT_TABULATED(ffreq,filter_xspect_avg))*100,prec=3)+' %',$
               'synth/cross power/40-400kHz: '+((fidx[0] LE 1)?'NaN':$
               rm0es(SQRT(INT_TABULATED(ffreq[fidx],filter_xspect_avg[fidx]))*100,prec=3)+' %')]

     FOR s = 0, N_ELEMENTS(infostr)-1 DO print, infostr[s]
     print, ''

     set_output, diag, sp, header=['f/kHz','raw autopow','raw crosspow.', $
                                   'synth autopow', 'synth crosspow'], $
                dat=[[ffreq],[REFORM(raw_aspect_avg)],[REFORM(raw_xspect_avg)],$
                     [REFORM(filter_aspect_avg)],[REFORM(filter_xspect_avg)]],$
                commentline=[[infostr],'df/kHz = '+rm0es(dffreq),$
                             'f_nyquist/kHz = '+rm0es(freq_nyquist)]
     ENDIF

     IF (radcorr) THEN BEGIN
        set_output, diag, sp, header=['lag/mm','raw corr','synth corr'], $
                    dat=[[lagr],[rawcorr],[filtercorr]],$
                    commentline=['raw FWHM [mm]   = '+rm0es(raw_corrfwhm),$
                                 'synth FWHM [mm] = '+rm0es(synth_corrfwhm),$
                                 'raw 1/e [mm]    = '+rm0es(raw_corrinve),$
                                 'synth 1/e [mm]  = '+rm0es(synth_corrinve)],$
                    append=freq_spect
     ENDIF

     set_output, diag, sp, /reset

  ENDFOR ;isp

END
