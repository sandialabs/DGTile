-- this is an example application input
-- 

input = {

  name = "advection",
  gamma = 1.4,

  time = {
    cfl = 0.9,
    end_time = 1.,
    to_terminal = {{kind='at_always'}}
  },

  basis = {
    polynomial_order = 1,
    quadrature_rule = 2,
    tensor_product = true
  },

  mesh = {
    X = {num_blocks = 2, num_cells = 10, min = 0., max = 1.},
    Y = {num_blocks = 2, num_cells = 10, min = 0., max = 1.},
  },

  initial_conditions = {
    density = function(x,y,z)
      xi = x + y + z
      return 1. + 0.1*math.sin(2*math.pi*xi)
    end,
    pressure = function(x,y,z)
      return 2.
    end,
    velocity = function(x,y,z)
      return 1.,1.,0.
    end
  },

}

dgtile.run(input)
