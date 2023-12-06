input = {

  name = 'advection',
  
  time = {
    cfl = 0.9,
    end_time = 1.,
    to_terminal = {{kind='at_step_periodically', frequency=20}},
    integrator = 'ssprk2'
  },

  basis = {
    polynomial_order = 1,
    quadrature_rule = 2,
    tensor_product = true
  },

  mesh = {
    X = {num_blocks = 2, num_cells = 20, min = 0., max = 1., periodic = true},
    Y = {num_blocks = 2, num_cells = 20, min = 0., max = 1., periodic = true}
  },

  hydro = {
    gamma = 1.4,
    density = function(x,y,z)
      xi = x + y + z
      return 1. + 0.1*math.sin(2*math.pi*xi)
    end,
    velocity = {1.,1.,1.},
    pressure = 2.
  },

  vtk = {
    when = {{kind='at_time_periodically', frequency=0.1}}
  }

}

dgtile.run(input)
