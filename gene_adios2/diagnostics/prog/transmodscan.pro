FUNCTION transmodscan_info

  COMMON global_vars

  IF PTR_VALID(scanlog) THEN BEGIN
    n_vars = N_ELEMENTS((*scanlog).dimarr)+4
    IF n_vars GT 0 THEN BEGIN
      ext_vars = STRARR(3,n_vars)
      ext_vars[0,0] = 'eigenvalues'
      ext_vars[0,n_vars-1] = 'exponent'
      ext_vars[0,n_vars-2] = 'kxind'
      ext_vars[0,n_vars-3] = 'const'
      ext_vars[0,n_vars-4] = 'nrgind'
      ext_vars[0,1:n_vars-5] = (*scanlog).dimarr[1:n_vars-5]
      ext_vars[1,*] = '0'
      ext_vars[2,0:n_vars-5] = 'index between '+rm0es(0)+' and '+$
        rm0es((*scanlog).dsize[*]-1)+' (corresponds to values from '+$
	rm0es(STRING([1,(*scanlog).axxis[*,0]]))+' to '+$
        rm0es(STRING([(*scanlog).dsize[0]-1,(*scanlog).rd[*]]))+'); default: * (all)'
      ext_vars[2,n_vars-1]='weighting function exponent; default 2'
      ext_vars[2,n_vars-2]='index of kx modes considered; default: -2 (all)'
      ext_vars[2,n_vars-3]='renomalization constant for transport values; default: 1.0'
      ext_vars[2,n_vars-4]='nrg indices for ql transport levels, first index sets reference; default: '+$
                           rm0es((par.n_spec-1)*par.nrgcols+6)+((par.n_spec GT 2)?',..,'+$
                           rm0es(par.nrgcols+6):'')+',6 (Q_es; last species is reference)'
    ENDIF ELSE ext_vars = ''
  ENDIF ELSE ext_vars = ''

  RETURN, {$
      type      : 'scan',$
      title     : 'Transport model',$
      help_text : ['Evaluates and plots the quasilinear transport model '+$
                   'const*(omt+omn)*gamma/<kperp^2> '+$
                   'for the first nrgind value. All remaining nrgind choices are derived from '+$
                   'this value via their linear ratios. In above model, <kperp^2> denotes '+$
                   'int k_perp^2 |phi|^exponent dx dz / int |phi|^exponent dx dz. '+$
                   'Sometimes such models just consider the ky component of k_perp which '+$
                   'is why it is additional written in the main output file together with '+$
                   'the growth rate gamma and all chosen scan parameters. '+$
                   'NOTE: Maximization or summation over ky still needs to be (re)implemented '+$
                   'and has to be performed by the user for now'],$
      ext_vars  : ext_vars}

END

;#####################################################################


PRO transmodscan_init,diag

  COMMON global_vars

  IF (*diag).table_entry[0,0] EQ '' THEN RETURN

  e = (*diag).external_vars
  IF NOT KEYWORD_SET(*e.exponent) THEN (*e.exponent)=2
  IF NOT KEYWORD_SET(*e.kxind) THEN (*e.kxind)=INDGEN(par.nkx0-1+(par.nkx0 MOD 2))
  IF NOT KEYWORD_SET(*e.const) THEN (*e.const)=1
  IF NOT KEYWORD_SET(*e.nrgind) THEN (*e.nrgind)=REVERSE(INDGEN(par.n_spec))*par.nrgcols+6
  n_vars = N_ELEMENTS((*diag).table_entry[3,*])
  IF n_vars GE 1 THEN eigenvalues = (*diag).table_entry[3,0]
  IF STRLEN(eigenvalues) NE 0 THEN BEGIN
     printerror, 'eigenvalues selection currently deactivated in transmodscan'
     eigenvalues='*'
  ENDIF
  IF STRLEN(eigenvalues) EQ 0 THEN eigenvalues='*'

  plotparA=REFORM((*diag).table_entry[3,1:n_vars-2])
  inds = WHERE(plotparA EQ '')
  IF (inds[0] GE 0) THEN plotparA[inds]='*'

  nevs = (*scanlog).dsize[0]
  nscans = PRODUCT((*scanlog).dsize[1:*],/INTEGER)
  ndsize = PRODUCT((*scanlog).dsize,/INTEGER)
  dummyarr = MAKE_ARRAY(nevs,nscans,/DOUBLE) ;(*scanlog).dsize,/DOUBLE)
  ratiofield=MAKE_ARRAY([N_ELEMENTS(*e.nrgind),nevs,nscans],/DOUBLE)

  maxsize=(*scanlog).dsize

  ind=1
  IF STRCMP(STRTRIM((*scanlog).dimarr[1],2),'ky',2) EQ 1 THEN maxsize(ind)=1

  max_ky=MAKE_ARRAY(maxsize,/DOUBLE)
  ky_max=MAKE_ARRAY(maxsize,/DOUBLE)

  o_nky0 = par.nky0-par.ky0_ind ;original nky0

  i = set_internal_vars(diag,{$
      eigenvalues  : eigenvalues,$
      nrgind       : *e.nrgind,$
      kxind        : *e.kxind,$
      const        : *e.const,$
      plotparA     : plotparA,$
      ky_max       : ky_max,$
      max_ky       : max_ky,$
      o_nky0       : o_nky0,$
      initrunnum   : 1,$
      gamofky      : DBLARR(o_nky0),$
      refsp        : (*e.nrgind)[0]/par.nrgcols,$
      ratiofield   : ratiofield,$
      gamma        : dummyarr,$
      qlkernel     : dummyarr,$ ;Kernel for QL model: gamma/<k_perp^2>
      qltrans      : dummyarr,$ ;Kernel*const*(omt_refsp+omn_refsp)
      kperp_sq     : dummyarr,$
      ky_sq        : dummyarr,$
      kyval        : dummyarr,$
      exponent     : *e.exponent,$
      max_arr      : DBLARR(o_nky0),$
      maxsize      : maxsize,$
      initrunlabel : series.run_labels})

  ; request field and nrg variables
  series.request_nrg = 1
  fft_format, kxky=0

END

;######################################################################
FUNCTION transmodscan_calcweight, diag,ix,iz
  COMMON global_vars

  i = (*diag).internal_vars

  iy = par.ky0_ind

  phi_exp = (ABS((*mom[0,0].kxky)[ix,iy,iz]))^(*i).exponent*$
            (*series.geom).jac_norm[iz]

  k_sq_weighed = (((*par.kx)[ix,iy])^2*$
              (*series.geom).gxx[iz]+(*par.ky)[iy]*(*par.kx)[ix,iy]*$
              2*(*series.geom).gxy[iz]+((*par.ky)[iy])^2*$
              (*series.geom).gyy[iz])*phi_exp

  ky_sq_weighed = ((*par.ky)[iy])^2*(*series.geom).gyy[iz]*phi_exp

  return, [k_sq_weighed, ky_sq_weighed, phi_exp]
END

;######################################################################
PRO transmodscan_calc, diag, inds, calc
  COMMON global_vars

  i = (*diag).internal_vars

  IF calc THEN BEGIN
      k_perp_sq = DBLARR((*i).o_nky0)
      ky_sq = DBLARR((*i).o_nky0)
      phi_sum = DBLARR((*i).o_nky0)
      FOR ikx = 0, N_ELEMENTS((*i).kxind)-1 DO BEGIN
         FOR iz = 0, par.nz0-1 DO BEGIN
            retval=transmodscan_calcweight(diag,(*i).kxind[ikx],iz)
            k_perp_sq=k_perp_sq + retval[0]
            ky_sq = ky_sq+retval[1]
            phi_sum=phi_sum +retval[2]
         ENDFOR
      ENDFOR
      k_perp_sq /= phi_sum
      ky_sq /= phi_sum

      (*i).gamma[inds] = REAL_PART((*scanlog).eigenfield[inds])
      (*i).qlkernel[inds]=REAL_PART((*scanlog).eigenfield[inds])/k_perp_sq
      (*i).qltrans[inds] = (*i).qlkernel[inds]*(*i).const*$
                           (spec[(*i).refsp].omt+spec[(*i).refsp].omn)
      (*i).kyval[inds] = (*par.ky)[par.ky0_ind:*]
      (*i).ky_sq[inds] = ky_sq
      (*i).kperp_sq[inds] = k_perp_sq

  ENDIF ELSE print, 'missing or stable eigenvalue'

END



;######################################################################

PRO transmodscan_loop, diag

  COMMON global_vars

  i = (*diag).internal_vars

  iscan=FIX(series.run_labels)
  nnrgind = N_ELEMENTS((*i).nrgind)

; check whether value is supposed to be considered for further post-processing
  calc=scan_plotnec(diag,[(*i).eigenvalues,(*i).plotparA],iscan)

  IF NOT calc THEN RETURN

; calc nrgindex for run and set ratiofield
  IF par.comp_type EQ 'IV'  THEN BEGIN
     ind_sortf = (iscan-1)*(*scanlog).dsize[0]
     IF (*scanlog).sortfield[ind_sortf] LT 0 THEN BEGIN
        ;select last nrg time step
        nrgtindex = N_ELEMENTS((*nrg).time) - 1
        FOR ind=0,nnrgind-1 DO BEGIN
           (*i).ratiofield[ind,0,iscan-1]=$
              (*nrg)[nrgtindex].data[(*i).nrgind[ind]]
        ENDFOR
        inds = INDGEN((*i).o_nky0)+iscan-1
        transmodscan_calc,diag,inds,calc
     ENDIF
  ENDIF ELSE BEGIN
     found=0
     FOR isort = 0, (*scanlog).dsize[0]-1 DO BEGIN
        ind_sortf = (iscan-1)*(*scanlog).dsize[0]+isort
        ind_ev = -(*scanlog).sortfield[ind_sortf]
        IF (ind_ev EQ FIX(mom_time)) THEN BEGIN
           ;select same nrg entry as in mom file
           nrgtindex = WHERE((*nrg).time EQ mom_time)
           FOR ind=0,nnrgind-1 DO BEGIN
              (*i).ratiofield[ind,ind_ev-1,iscan-1]=$
                 (*nrg)[nrgtindex].data[(*i).nrgind[ind]]
           ENDFOR
           transmodscan_calc,diag,ind_sortf,calc
           found+=1
        ENDIF
     ENDFOR
     IF found LT 1 THEN print,'stable or no EVs found'
  ENDELSE

IF (1 EQ 0) THEN BEGIN
  IF par.comp_type EQ 'IV' THEN BEGIN ; initvalues
;    if scan over ky --> max after scanrange
     IF STRCMP(STRTRIM((*scanlog).dimarr[1],2),'ky',2) EQ 1 THEN BEGIN
        IF (FIX(series.run_labels) MOD ((*scanlog).dsize[1]/(*i).o_nky0)$
            EQ 0) AND (FIX(series.run_labels) GT 0) THEN BEGIN
           (*i).max_ky[(*i).initrunnum-1]=MAX((*i).$
           qlkernel[((*i).initrunnum-1)*(*i).o_nky0*$
                    (*scanlog).dsize[1]:(*i).initrunnum*(*i).o_nky0*$
                    (*scanlog).dsize[1]-1],Max_Subscript,/NAN)
           (*i).ky_max[(*i).initrunnum-1]=(*scanlog).rd[0]-$
                     ((*scanlog).dsize[1]-Max_Subscript-1)*(*scanlog).deltad[0]
           (*i).initrunnum+=1
        ENDIF
;    max of ky if nky > 1
     ENDIF ELSE BEGIN
        (*i).max_ky[(*i).initrunnum-1]=$
           MAX((*i).qlkernel[(*i).initrunnum*(*i).o_nky0-1:$
           (*i).initrunnum*(*i).o_nky0+(*i).o_nky0-2],Max_Subscript,/NAN)
        (*i).ky_max[(*i).initrunnum-1]=(*scanlog).rd[0]-$
           ((*scanlog).dsize[1]-Max_Subscript-1)*(*scanlog).deltad[0]
        (*i).initrunnum+=1
     ENDELSE
  ENDIF ELSE BEGIN ; eigenvalues
;    if last mom entry in this particular run
     IF FIX(mom_time) EQ (*scanlog).dsize[0] THEN BEGIN
;       if scanned over ky
        IF STRCMP(STRTRIM((*scanlog).dimarr[1],2),'ky',2) EQ 1 THEN BEGIN
;          if end of ky-direction (and not the first run)
           IF iscan MOD (*scanlog).dsize[1] EQ 0 $
              AND iscan GT 1 THEN BEGIN
;             max of ky-range
              FOR iev=0, (*scanlog).dsize[0]-1 DO BEGIN
                 (*i).max_ky[(iscan/(*scanlog).$
                    dsize[1]-1)*(*scanlog).dsize[0]+iev]=$
                    MAX((*i).qlkernel[((iscan-$
                    (*scanlog).dsize[1])*(*scanlog).dsize[0]+iev):$
                    iscan*(*scanlog).dsize[0]-1:$
                    ((*scanlog).dsize[0])],Max_Subscript,/NAN)
                 (*i).ky_max[iscan/(*scanlog).$
                    dsize[1]-1]=(*scanlog).rd[0]-$
                    ((*scanlog).dsize[1]-Max_Subscript-1)*$
                    (*scanlog).deltad[0]
                 ENDFOR
           ENDIF
        ENDIF ELSE BEGIN
           FOR iev=0, (*scanlog).dsize[0]-1 DO BEGIN
              ind_ev = (iscan-1)*(*scanlog).dsize[0]+iev
              (*i).max_ky[ind_ev]=(*i).qlkernel[ind_ev]
           ENDFOR
        ENDELSE
     ENDIF
  ENDELSE
ENDIF

END



;######################################################################

PRO transmodscan_output, diag

  COMMON global_vars

  i = (*diag).internal_vars
  series.run_labels = (*i).initrunlabel

  nnrgind = N_ELEMENTS((*i).nrgind)
  ndsize = PRODUCT((*scanlog).dsize,/INTEGER)
  ndims = N_ELEMENTS((*scanlog).dimarr)

  kyminind = WHERE((*scanlog).dimarr EQ 'kymin')

;*************calculate quasilinear ratios***************************
  ratiofield = REFORM((*i).ratiofield,nnrgind,ndsize)
  FOR ind=nnrgind-1,0,-1 DO $
      ratiofield[ind,*]=ratiofield[ind,*]/ratiofield[0,*]*(*i).qltrans
  (*i).ratiofield = ratiofield ;REFORM(ratiofield,[nnrgind,[(*scanlog).dsize]],/OVERWRITE)

;********************************************************************

  iev = 0
  kind = WHERE((*i).kyval[iev,*] GT 0)

  header = ['ky','<ky^2 gyy>','<kperp^2>']
  data = [[REFORM((*i).kyval[iev,kind])],[REFORM((*i).ky_sq[iev,kind])],$
          [REFORM((*i).kperp_sq[iev,kind])]]
  labelstr = ''
  ichar = 120B
  FOR s = 1, ndims-1 DO BEGIN
     IF (s NE kyminind) THEN BEGIN
        header = [header,(*scanlog).dimarr[s]]
        data = [[data],[REFORM((*scanlog).scanpars[s-1,kind])]]
     ENDIF
     IF (*i).plotparA[s-1] EQ '*' THEN BEGIN
        labelstr += 'set '+STRING(ichar)+'label "'+(*scanlog).dimarr[s]+'"; '
        ichar += 1B
     ENDIF
  ENDFOR

  header = [header,'gamma','qlkernel']
  data = [[data],[REFORM((*i).gamma[iev,kind])],[REFORM((*i).qlkernel[iev,kind])]]

  FOR n = 0, nnrgind-1 DO BEGIN
     header = [header,get_nrg_string((*i).nrgind[n],/addspec)]
     data = [[data],[REFORM((*i).ratiofield[n,iev,kind])]]
  ENDFOR

  set_output, diag, header = header,$
              commentlines=['reference species: '+spec[(*i).refsp].name,$
                            'ky : original ky values where gamma>0',$
                           '<ky^2 gyy> = int ky^2 gyy |phi|^'+rm0es((*i).exponent)+$
                            ' J dx dz/int |phi|^'+rm0es((*i).exponent)+' J dx dz',$
                           '<kperp^2> = int (gxx kx^2 + gxy kx ky + gyy ky^2) |phi|^'+$
                            rm0es((*i).exponent)+' J dx dz/int |phi|^'+$
                            rm0es((*i).exponent)+' J dx dz',$
                           'gamma : growth rate from IV or EV run',$
                           'qlkernel : quasilinear kernel gamma/<kperp^2>',$
                           get_nrg_string((*i).nrgind[0],/addspec)+$
                            ' : const*(omt+omn)*qlkernel','',$
                           (par.comp_type EQ 'EV')?'EV 1 (POSSIBLY REARRANGED)':''],$
                           dat=data

  FOR iev = 1, (*scanlog).dsize[0]-1 DO BEGIN
     kind = WHERE(REFORM((*i).kyval[iev,*]) GT 0)
     IF kind[0] NE -1 THEN BEGIN
        data = [[REFORM((*i).kyval[iev,kind])],[REFORM((*i).ky_sq[iev,kind])],$
                [REFORM((*i).kperp_sq[iev,kind])]]
        FOR s = 1, ndims-1 DO IF (s NE kyminind) THEN $
              data = [[data],[REFORM((*scanlog).scanpars[s-1,kind])]]
        data = [[data],[REFORM((*i).gamma[iev,kind])],[REFORM((*i).qlkernel[iev,kind])]]
        FOR n = 0, nnrgind-1 DO data = [[data],[REFORM((*i).ratiofield[n,iev,kind])]]
        set_output, diag, header = header,$
              commentlines=[(par.comp_type EQ 'EV')?'EV '+rm0es(iev+1)+' (POSSIBLY REARRANGED)':''],$
              dat=data,/APPEND
     ENDIF
  ENDFOR

;************************plot all ratio fields***********************

  FOR n=0,nnrgind-1 DO BEGIN
     scan_plot, diag, 'quasilinear transport',REFORM(ratiofield[n,*]),$
                'abs', 0, 1, 2, 0.01, 0.06, 0.99, 0.16, 2,$
                suffix='-rval_'+rm0es(STRING(n))

     IF gui.out.out_format[0] AND ichar LT 123B THEN SPAWN,"echo "+"'"+labelstr+$
               'set '+STRING(ichar)+'label "'+$
               get_nrg_string((*i).nrgind[n],/addspec)+'"; splot '+$
               '"'+gui.out.out_path+'transmodscan-rval_'+rm0es(STRING(n))+$
               '_0001.dat"'+" u 1:2:3 w lp' | gnuplot -persist"
  ENDFOR

RETURN

;************************plot maximization over ky*********************
  IF NOT ((*scanlog).dsize[1] EQ (*i).maxsize[1]) THEN BEGIN
      IF STRCMP(STRTRIM((*scanlog).dimarr[1],2),'ky',2) EQ 1 $
        THEN plotparA1='0'
      dummy=(*scanlog).dsize
      (*scanlog).dsize=(*i).maxsize

      maxplotfield=(*i).max_ky
      ;printjob for all given nrgind
      FOR ind=0,N_ELEMENTS((*i).nrgind)-1 DO BEGIN
          FOR subs=0,N_ELEMENTS((*i).max_ky)-1 DO BEGIN
              adummy=(*i).max_ky[subs]*$
                     (*i).ratiofield[[ind+N_ELEMENTS((*i).nrgind)*WHERE((*i).qlkernel EQ (*i).max_ky[subs])]]
              maxplotfield[subs]=adummy[0]
          ENDFOR
          scan_plot, diag, 'max of quasilinear transport',maxplotfield,'abs',$
                     0, 1, 3, 0.01, 0.06, 0.99, 0.16, 2, suffix='m-rval_'+rm0es(STRING(ind))
          IF gui.out.out_format[0] THEN SPAWN,"echo "+"'"+'set ylabel "nrg '+rm0es(STRING((*i).nrgind[ind]))+$
                '"; plot '+'"'+gui.out.out_path+'transmodscanm-rval_'+$
                rm0es(STRING(ind))+'_0001.dat"'+" w lp' | gnuplot -persist"
      ENDFOR
      scan_plot, diag, 'ky_max',(*i).ky_max,'abs',$
                 0, 1, 3, 0.01, 0.06, 0.99, 0.16, 2, suffix='kymax'
      (*scanlog).dsize=dummy
  ENDIF



END
