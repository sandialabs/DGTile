whens = {
  [1] = {{kind='at_always'}},
  [2] = {{kind='at_never'}},
  [3] = {{kind='at_step', step=3}},
  [4] = {{kind='at_step', step=2},{kind='at_step', step=4}},
  [5] = {
    {kind='at_step', step=0},
    {kind='at_step', step=1},
    {kind='at_step', step=2},
    {kind='at_step', step=3},
    {kind='at_step', step=4},
    {kind='at_step', step=5}
  },
  [6] = {{kind='at_step_periodically', frequency=3}},
  [7] = {{kind='at_time', time=0.15}},
  [8] = {{kind='at_exact_time', time=0.15}},
  [9] = {{kind='at_time_periodically', frequency=0.2}},
  [10] = {{kind='at_exact_time_periodically', frequency=0.15}}
}

dgt_check_whens(whens)
