
#RAPT algorithm parameters from Talkin, D., & Kleijn, W. B. (1995). A robust algorithm for pitch tracking (RAPT). Speech coding and synthesis, 495, 518.
const F0min = 50.0
const F0max = 500.0
const frame_step = 0.01
const nccf_win_size = 0.0075
const CAND_TR = 0.3
const LAG_WT = 0.3
const FREQ_WT = 0.02
const VTRAN_C = 0.005
const VTR_A_C = 0.5
const VTR_S_C = 0.5
const VO_BIAS = 0.0
const DOUBL_C = 0.35
const A_FACT = 10000.0
const N_CANDS = 20

