# ifndef _WAVE_H_
# define _WAVE_H_

# define  P_Wave   0
# define  SV_Wave  1
# define  SH_Wave  2

typedef struct {
  int   type;
  double slowness;
} Wave;

# endif
