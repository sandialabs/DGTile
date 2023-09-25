-- this is an example application input
-- 

input = {

  name = "advection",
  num_materials = 1,

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

  materials = {
    gamma_0 = 1.4
  },

  initial_conditions = {
    density_0 = function(x,y,z)
      xi = x + y + z
      return 1. + 0.1*math.sin(2*math.pi*xi)
    end,
    pressure_0 = function(x,y,z)
      return 2.
    end,
    velocity = function(x,y,z)
      return 1.,1.,1.
    end
  },

}

dgtile.run(input)
