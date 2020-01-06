FUNCTION nrgsummary_info

  RETURN, {$
    type      : 'nrg',$
    title     : 'summary',$
    help_text : ['Plots all nrg data, lists growth rates and'+$
                 'saturation levels'],$
    ext_vars  : [['t_corr','0','correlation time for error '+$
                  'computation; default: -1 (automatic calculation) '+$
                  'for nonlinear, 0 (no errors) for linear runs; '+$
                  'set to positive value for manual setting'],$
                 ['show_SI','1','include table with SI values'],$
                 ['show_sdev','1','include standard deviation '+$
                  'with saturated values']]}

END

;######################################################################

PRO nrgsummary_init, diag

END

;######################################################################

PRO nrgsummary_loop, diag

END

;######################################################################

PRO nrgsummary_output, diag

  COMMON global_vars

  e = (*diag).external_vars
  t_corr = N_ELEMENTS(*e.t_corr) NE 1 ? -par.nonlinear : *e.t_corr
  show_SI = KEYWORD_SET(*e.show_SI)
  show_sdev = KEYWORD_SET(*e.show_sdev)

  plot_nrg_data, diag=diag, t_corr_in=t_corr, show_SI=show_SI, $
    show_sdev=show_sdev
  
END
