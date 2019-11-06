
#ifdef __cplusplus
extern "C" {
#endif

  int ol_register_process(const char *process, int amptype);
  void ol_start();
  void ol_finish();

  void ol_setparameter_int(const char *key, int val);
  void ol_getparameter_int(const char *key, int *val);
  void ol_setparameter_double(const char *key, double val);
  void ol_getparameter_double(const char *key, double *val);
  void ol_setparameter_string(const char *key, const char *val);

  void ol_evaluate_tree(int id, const double *pp, double *m2l0);
  void ol_tree_colbasis_dim(int id, int *ncolb, int *colelemsz, int *nhel);
  void ol_tree_colbasis(int id, int *basis, int *needed);
  void ol_tree_colourflow(int id, int *flowbasis);
  void ol_evaluate_tree_colvect(int id, const double *pp, double *amp, int *nhel);
  void ol_evaluate_tree_colvect2(int id, const double *pp, double *m2arr);
  void ol_evaluate_associated(int id, const double *pp, int level, double *res);

  void ol_evaluate_full(int id, const double *pp, double *m2l0, double *m2l1, double *ir1, double *m2l2, double *ir2, double *acc);
  void ol_evaluate_loop(int id, const double *pp, double *m2l0, double *m2l1, double *acc);
  void ol_evaluate_loop2(int id, const double *pp, double *m2l0, double *acc);
  void ol_evaluate_sc(int id, const double *pp, int emitter, double *polvect, double *m2sc);
  void ol_evaluate_sc2(int id, const double *pp, int emitter, double *polvect, double *m2sc);
  void ol_evaluate_loopsc(int id, const double *pp, int emitter, double *polvect, double *m2sc);
  void ol_evaluate_cc(int id, const double *pp, double *tree, double *m2cc, double *m2ewcc);
  void ol_evaluate_cc2(int id, const double *pp, double *tree, double *m2cc, double *m2ewcc);
  void ol_evaluate_ccmatrix(int id, const double *pp, double *tree, double *m2cc, double *m2ewcc);
  void ol_evaluate_ccmatrix2(int id, const double *pp, double *tree, double *m2cc, double *m2ewcc);
  void ol_evaluate_ccewmatrix(int id, const double *pp, double *tree, double *m2ccew);
  void ol_evaluate_ccewmatrix2(int id, const double *pp, double *tree, double *m2ccew);
  void ol_evaluate_loopcc(int id, const double *pp, double *tree, double *m2l1, double *cc, double *ewcc);
  void ol_evaluate_loopccmatrix(int id, const double *pp, double *tree, double *m2l1, double *ccij, double *ewcc);
  void ol_evaluate_scpowheg(int id, const double *pp, int emitter, double *res, double *resmunu);
  void ol_evaluate_sctensor(int id, const double *pp, int emitter, double *res, double *resmunu);
  void ol_evaluate_scpowheg2(int id, const double *pp, int emitter, double *res, double *resmunu);
  void ol_evaluate_sctensor2(int id, const double *pp, int emitter, double *res, double *resmunu);
  void ol_evaluate_loopscpowheg(int id, const double *pp, int emitter, double *res, double *resmunu);
  void ol_evaluate_loopsctensor(int id, const double *pp, int emitter, double *res, double *resmunu);
  void ol_evaluate_iop(int id, const double *pp, double *res, double *resir, double *acc);
  void ol_evaluate_iop2(int id, const double *pp, double *res, double *resir, double *acc);
  void ol_evaluate_ct(int id, const double *pp, double *m2l0, double *m2ct);
  void ol_evaluate_loopct(int id, const double *pp, double *m2l0, double *m2ct);
  void ol_evaluate_r2(int id, const double *pp, double *m2l0, double *m2r2);
  void ol_evaluate_pt(int id, const double *pp, double *m2l0, double *m2pt, double *m2l1);
  void ol_evaluate_poles(int id, const double *pp, double *m2l0, double *m2bare, double *m2ct, double *m2ir, double *m2sum);
  void ol_evaluate_m2schsf(int id, const double *pp, double *m2l0, double *m2schsf);

  void ol_welcome(char *str);
  void ol_version_string(char *str);
  void ol_set_init_error_fatal(int flag);
  int  ol_get_error();
  int ol_n_external(int id);
  void ol_phase_space_point(int id, double sqrt_s, double *pp);
  void ol_printparameter(char *filename);
  void ol_parameters_flush();
  void ol_tree_parameters_flush();
  int ol_amplitudetype(int id);

  // BLHA interface
  void OLP_SetParameter(const char *para, const double *re, const double *im, int *ierr);
  void OLP_EvalSubProcess(const int *id, const double *pp, const double *mu, const double *alpha_s, double *rval);
  void OLP_EvalSubProcess2(const int *id, const double *pp, const double *mu, double *rval, double *acc);
  void olp_scpolvec(const int *emitter, const int *mom);
  void OLP_Info(char *olp_name, char *olp_version, char *message);
  void OLP_PrintParameter(const char *filename);
  void OLP_Start(const char *contract_file_name, int *ierr);
  void OLP_StartLine(const char *contract_line, char *answer_line, int *ierr);

#ifdef __cplusplus
}
#endif
