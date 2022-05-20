*********************************************************************
*                PYQM CONTROL CARD Updated 11-05-2022               *
*********************************************************************
*                                                                   *
*               control card:  codewd = PYQM-CTRL                   *
*                                                                   *
*       Added by Mark to give more options for PYQM                 *
*                                                                   *
*       what (1) = fraction of PYQM recoil applied to nucleus (E*)  *
*                                                   default=0.0     *
*       what (2) = PyQM Pt model (iPtF in ApplyQW)                  *
*                  0: No Pt broadening                              *
*                  1: DPt2 = qhat * length                          *
*                  2: DPt2 proportional to energy loss QW_w         *
*                  3: No Pt broadening due to collinear gluons      *
*                                                                   *
*       what (3) = Emit recoil as a single gluon?(gluons)           *
*                  0: No gluons (def)                               *
*                  1: single hard gluon                             *
*                  2: one hard gluon + soft                         *
*                  3: soft gluons                                   *
*                                                                   *
*       what (4) = SupFactor (needed for iPtF=2)  default=1.0       *
*                                                                   *
*       what (5) = New SW calculation for heavy quarks(default=1.0) *
*                  0: No                                            *  
*                  1: Yes                                           *
*                                                                   *
*       what (6) = Energy parton threshold (iet, default=0.25 GeV)  * 
*      (This line control should be manipulated only  by developers)*
*                                                                   *
*********************************************************************
Example:
To enable the PYQM control card, go to the input file muXe_pyqm_qhat_g1.inp and 
manipulate the individual options, although we do not recommend that you change 
"iet" as a user.

Option 1, generating a single hard gluon calculating DPt2 as qhat * length:
 
*PYQM_CTRL E*  Pt model #gluons Supfac SW    iet  
PYQM-CTRL  0.0    1.0    1.0    1.0    1.0   0.25


qhat can be manipulated in the Pythia control card, such as S1ALL003,which is located
under Examples folder in the main directory.

